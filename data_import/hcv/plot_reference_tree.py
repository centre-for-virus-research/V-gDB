from __future__ import annotations

import argparse
import csv
import html
import json
import sqlite3
import sys
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_highlight_tree import PALETTE, TipIndex, read_metadata  # noqa: E402

HERE = Path(__file__).resolve().parent
OUTGROUPS = ["NC_009528.2", "NC_009527.1", "NC_020807.1"]
# Same major-clade colours as Analysis/Plot_trees; Africa-3 gets the spare slot.
MAJOR_COLOR = {
    "Cosmopolitan": "#2a78d6", "Bats": "#eb6834", "Asian": "#1baf7a", "Arctic": "#eda100",
    "RAC-SK": "#e87ba4", "Africa-2": "#008300", "Africa-3": "#4a3aa7", "Indian-Sub": "#e34948",
}
UNKNOWN = "#8a8984"
UNCLASSIFIED = "Unclassified"
BANNER = "MADDOG-V2 workflow deployment to the server is in testing mode"
META_COLUMNS = ["country", "host", "collection_year"]
EMPTY = {"", "null", "none", "nan"}


# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value.lower() in EMPTY else value


def read_two_columns(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        return {row[0].strip(): row[1].strip() for row in reader if len(row) >= 2}


def tip_names(path):
    """Tip names from a Newick file, or one name per line from a text file."""
    text = Path(path).read_text(encoding="utf-8")
    if text.lstrip().startswith("("):
        return {t.name for t in Phylo.read(path, "newick").get_terminals() if t.name}
    return {line.split()[0] for line in text.splitlines() if line.strip() and not line.startswith("#")}


def read_clades(tips, major_path, minor_path, db):
    """{tip: (major, minor)}: EPA tables of the references first, then the database."""
    out = {}
    if major_path and minor_path and Path(major_path).exists() and Path(minor_path).exists():
        major, minor = read_two_columns(major_path), read_two_columns(minor_path)
        for acc in major:
            tip = tips.find(acc)
            if tip:
                out[tip] = (clean(major[acc]), clean(minor.get(acc)))
    if db:
        con = sqlite3.connect(db)
        for acc, version, ma, mi in con.execute(
                "SELECT primary_accession, accession_version, EPA_major_clade, EPA_minor_clade FROM meta_data"):
            tip = tips.find(acc or "") or tips.find(version or "")
            if tip and tip not in out and clean(ma):
                out[tip] = (clean(ma), clean(mi))
        con.close()
    return out


def read_lineages(lineage_dir, tip_sub, tips):
    """{tip: (lineage, parent)} from each sub-clade's lineage table; a sequence only takes
    the lineage from the table of its own sub-clade."""
    out, missing = {}, []
    for sub in sorted(set(tip_sub.values())):
        if sub == UNCLASSIFIED:
            continue
        base = Path(lineage_dir) / sub
        path = next((p for p in (base / "meta_data.tsv", base / f"{sub}.tip_lineage_assignments.tsv")
                     if p.exists()), None)
        if path is None:
            missing.append(sub)
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                tip = tips.find(row.get("ID") or row.get("sequence_id") or "")
                if tip and tip_sub.get(tip) == sub:
                    out[tip] = (clean(row.get("Lineage") or row.get("lineage")),
                                clean(row.get("Parent") or row.get("lineage_parent")))
    if missing:
        print(f"Warning: no lineage table in {lineage_dir} for: {', '.join(missing)}")
    return out


def read_db_metadata(db, tips):
    con = sqlite3.connect(db)
    values = {c: {} for c in META_COLUMNS}
    query = "SELECT primary_accession, accession_version, {} FROM meta_data".format(", ".join(META_COLUMNS))
    for acc, version, *row in con.execute(query):
        tip = tips.find(acc or "") or tips.find(version or "")
        if tip is None:
            continue
        for column, value in zip(META_COLUMNS, row):
            value = clean(value)
            if value:
                values[column].setdefault(tip, value)
    con.close()
    return values


# ------------------------------------------------------------
# Tree handling
# ------------------------------------------------------------

def subset_tree(tree, keep):
    """Tree reduced to the tips in `keep`, root kept, branch lengths summed through
    removed nodes. Iterative, so large trees are fine."""
    Clade = Phylo.BaseTree.Clade
    order, stack = [], [tree.root]
    while stack:
        clade = stack.pop()
        order.append(clade)
        stack.extend(clade.clades)
    new = {}
    for clade in reversed(order):
        if clade.is_terminal():
            if clade.name in keep:
                new[clade] = Clade(branch_length=clade.branch_length, name=clade.name)
            continue
        kids = [new[c] for c in clade.clades if c in new]
        if not kids:
            continue
        if len(kids) == 1 and clade is not tree.root:
            kid = kids[0]
            kid.branch_length = (kid.branch_length or 0.0) + (clade.branch_length or 0.0)
            new[clade] = kid
            continue
        node = Clade(branch_length=clade.branch_length, name=clade.name, clades=kids)
        node.confidence = clade.confidence
        new[clade] = node
    root = new.get(tree.root)
    if root is None:
        raise ValueError("none of the requested tips are in the tree")
    while len(root.clades) == 1:
        root = root.clades[0]
    root.branch_length = 0.0
    return Phylo.BaseTree.Tree(root=root, rooted=True)


def root_tree(tree, outgroups, keep):
    names = {t.name for t in tree.get_terminals()}
    present = [o for o in outgroups if o in names]
    if not present:
        print("No outgroup in the tree; keeping the root from the file")
        return tree
    tree.root_with_outgroup(*present)
    if not keep:
        tree = subset_tree(tree, names - set(present))
    return tree


def flat_nodes(tree):
    """Ladderized preorder list of {"c", "b", "n"/"s"}; a subtree is a contiguous range."""
    order, stack = [], [tree.root]
    while stack:
        clade = stack.pop()
        order.append(clade)
        stack.extend(clade.clades)
    count = {}
    for clade in reversed(order):
        count[clade] = 1 if clade.is_terminal() else sum(count[c] for c in clade.clades)
    for clade in order:
        clade.clades.sort(key=lambda c: count[c])

    nodes, stack = [], [(tree.root, -1)]
    while stack:
        clade, parent = stack.pop()
        i = len(nodes)
        node = {"c": [], "b": float(f"{clade.branch_length or 0.0:.6g}")}
        if clade.is_terminal():
            node["n"] = clade.name or ""
        else:
            support = clade.confidence if clade.confidence is not None else clade.name
            if support is not None:
                node["s"] = str(support)
        nodes.append(node)
        if parent >= 0:
            nodes[parent]["c"].append(i)
        for child in reversed(clade.clades):
            stack.append((child, i))
    for node in reversed(nodes):
        node["nt"] = 1 if not node["c"] else sum(nodes[c]["nt"] for c in node["c"])
    end = [0] * len(nodes)
    size = [1] * len(nodes)
    for i in range(len(nodes) - 1, -1, -1):
        size[i] += sum(size[c] for c in nodes[i]["c"])
        end[i] = i + size[i]
    return nodes, end


def dominant(nodes, start, end, label):
    """{node: (most common tip label, its count)} for the preorder range [start, end)."""
    out, pending = {}, {}
    for i in range(end - 1, start - 1, -1):
        kids = nodes[i]["c"]
        if not kids:
            counts = {label(i): 1}
        else:
            parts = sorted((pending.pop(c) for c in kids), key=len, reverse=True)
            counts = parts[0]
            for other in parts[1:]:
                for key, value in other.items():
                    counts[key] = counts.get(key, 0) + value
        pending[i] = counts
        out[i] = max(counts.items(), key=lambda kv: kv[1])
    return out


def find_points(nodes, root, dom, purity):
    """Topmost internal nodes whose most common label covers >= purity of their tips."""
    points, stack = [], [root]
    while stack:
        i = stack.pop()
        if not nodes[i]["c"]:
            continue
        label, count = dom[i]
        if label is not None and count >= purity * nodes[i]["nt"]:
            points.append(i)
        else:
            stack.extend(nodes[i]["c"])
    return points


def lineage_prefix(lineage, depth):
    if not lineage:
        return None
    if lineage == "Unassigned":
        return lineage
    return ".".join(lineage.split(".")[:depth])


def annotate(nodes, end, tip_major, tip_sub, lineages, purity):
    """Per node: dominant major (mj) and sub-clade (sc), and the collapse points:
    pm = major clade, pt = sub-clade or lineage group (ps sub-clade, lin lineage, dp %)."""
    majors = sorted({m for m in tip_major.values()}, key=lambda m: (m not in MAJOR_COLOR, m == UNCLASSIFIED, m))
    subs = sorted(set(tip_sub.values()))
    major_index = {m: i for i, m in enumerate(majors)}
    sub_index = {s: i for i, s in enumerate(subs)}
    tips = [i for i, n in enumerate(nodes) if not n["c"]]
    for i in tips:
        name = nodes[i]["n"]
        nodes[i]["mj"] = major_index[tip_major[name]]
        nodes[i]["sc"] = sub_index[tip_sub[name]]

    n = len(nodes)
    dom_major = dominant(nodes, 0, n, lambda i: nodes[i]["mj"])
    dom_sub = dominant(nodes, 0, n, lambda i: nodes[i]["sc"])
    for i, node in enumerate(nodes):
        if node["c"]:
            for key, dom in (("mj", dom_major), ("sc", dom_sub)):
                label, count = dom[i]
                node[key] = label if count >= purity * node["nt"] else -1

    def mark(i, count):
        nodes[i]["pt"] = 1
        share = round(100 * count / nodes[i]["nt"])
        if share < 100:
            nodes[i]["dp"] = share
        else:
            nodes[i].pop("dp", None)

    for i in find_points(nodes, 0, dom_major, purity):
        nodes[i]["pm"] = 1

    max_depth = max((len(v[0].split(".")) for v in lineages.values() if v[0]), default=0)
    tip_lineage = {i: lineages.get(nodes[i]["n"], ("", ""))[0] for i in tips}

    def lineage_points(region, sub, depth, parent_label):
        if depth > max_depth or nodes[region]["nt"] < 3:
            return
        dom = dominant(nodes, region, end[region],
                       lambda i: lineage_prefix(tip_lineage.get(i), depth))
        for q in find_points(nodes, region, dom, purity):
            label, count = dom[q]
            if label != parent_label:
                mark(q, count)
                nodes[q]["ps"], nodes[q]["lin"] = sub, label
            lineage_points(q, sub, depth + 1, label)

    sub_points = find_points(nodes, 0, dom_sub, purity)
    for i in sub_points:
        mark(i, dom_sub[i][1])
        nodes[i]["ps"] = dom_sub[i][0]
        lineage_points(i, dom_sub[i][0], 1, None)
    nodes[0].pop("pt", None)
    nodes[0].pop("pm", None)

    sub_rows = []
    for s in subs:
        members = [i for i in tips if nodes[i]["sc"] == sub_index[s]]
        sub_rows.append({"name": s, "major": nodes[members[0]]["mj"], "n": len(members)})
    major_rows = [{"name": m, "color": MAJOR_COLOR.get(m, UNKNOWN if m == UNCLASSIFIED else PALETTE[(i + 8) % len(PALETTE)]),
                   "n": sum(1 for v in tip_major.values() if v == m)} for i, m in enumerate(majors)]
    print(f"Collapse groups: {len(sub_points)} sub-clade, "
          f"{sum(1 for x in nodes if x.get('pt')) - len(sub_points)} lineage (purity {purity})")
    return sub_rows, major_rows


def encode_tips(nodes, lineages, meta_values):
    """Lineage and metadata as indexes into shared value lists (keeps the page small)."""
    lin_list, lin_index = [""], {"": 0}
    meta_lists = [[""] for _ in META_COLUMNS]
    meta_index = [{"": 0} for _ in META_COLUMNS]

    def code(value, table, index):
        if value not in index:
            index[value] = len(table)
            table.append(value)
        return index[value]

    for node in nodes:
        if node["c"]:
            continue
        name = node["n"]
        lineage, parent = lineages.get(name, ("", ""))
        if lineage:
            node["li"] = code(lineage, lin_list, lin_index)
        if parent:
            node["lp"] = code(parent, lin_list, lin_index)
        m = [code(meta_values.get(c, {}).get(name, ""), meta_lists[f], meta_index[f])
             for f, c in enumerate(META_COLUMNS)]
        if any(m):
            node["m"] = m
    return lin_list, meta_lists


def write_html(output, data):
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    page = HTML_TEMPLATE.replace("__TITLE__", html.escape(data["title"])).replace("__DATA__", payload)
    output = Path(output)
    if output.suffix.lower() != ".html":
        output = output.with_suffix(".html")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding="utf-8")
    return output


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  :root {
    --surface: #fcfcfb; --text: #0b0b0b; --text-2: #52514e; --muted: #8a8984;
    --border: #e4e3df; --hl: #c2185b;
    --banner-bg: #fff4d6; --banner-text: #7a4b00; --banner-border: #f0c96a;
  }
  * { box-sizing: border-box; }
  html, body { margin: 0; height: 100%; }
  body {
    background: var(--surface); color: var(--text);
    font-family: "DejaVu Sans", Verdana, "Segoe UI", sans-serif; font-size: 13px;
    display: flex; flex-direction: column;
  }
  .banner {
    flex: none; padding: 7px 16px; text-align: center; font-weight: 600; font-size: 13px;
    background: var(--banner-bg); color: var(--banner-text); border-bottom: 1px solid var(--banner-border);
  }
  header { flex: none; padding: 10px 16px 8px; border-bottom: 1px solid var(--border); max-height: 52vh; overflow: auto; }
  h1 { font-size: 17px; font-weight: 600; margin: 0 0 2px; }
  .stats { color: var(--text-2); font-size: 12px; margin-bottom: 8px; }
  .finders { display: flex; flex-wrap: wrap; gap: 8px 22px; align-items: flex-end; margin-bottom: 8px; }
  .finder { display: flex; flex-direction: column; gap: 3px; }
  .finder > span { font-size: 11px; font-weight: 600; color: var(--text-2); text-transform: uppercase; letter-spacing: .03em; }
  .finder .row { display: flex; gap: 6px; align-items: center; }
  input[type=search] {
    font: inherit; padding: 5px 9px; width: 250px; max-width: 70vw;
    border: 1px solid #cfcec9; border-radius: 5px; background: #fff; color: var(--text);
  }
  input[type=search]:focus { outline: 2px solid #2a78d6; outline-offset: -1px; }
  .msg { color: var(--muted); font-size: 12px; }
  .msg.err { color: #b3261e; }
  .controls { display: flex; flex-wrap: wrap; gap: 6px 14px; align-items: center; }
  .controls label { display: inline-flex; align-items: center; gap: 5px; color: var(--text-2); font-size: 12px; }
  button {
    font: inherit; font-size: 12px; padding: 4px 9px; cursor: pointer;
    border: 1px solid #cfcec9; border-radius: 5px; background: #fff; color: var(--text);
  }
  button:hover { background: #f1f0ed; }
  input[type=range] { width: 90px; }
  #card {
    display: none; margin-top: 8px; padding: 8px 12px; border: 1px solid var(--border);
    border-left: 4px solid var(--hl); border-radius: 6px; background: #fff;
  }
  #card .k { color: var(--muted); font-size: 11px; text-transform: uppercase; letter-spacing: .03em; }
  #card .v { font-weight: 600; font-size: 13px; }
  #card .grid { display: flex; flex-wrap: wrap; gap: 6px 22px; }
  #card .close { float: right; border: none; background: none; color: var(--muted); font-size: 16px; padding: 0 2px; }
  details.legend { margin-top: 8px; }
  details.legend summary { cursor: pointer; font-size: 12px; font-weight: 600; }
  .legend-items { display: flex; flex-wrap: wrap; gap: 4px 16px; margin-top: 6px; }
  .legend-item { display: inline-flex; align-items: center; gap: 6px; font-size: 12px; cursor: pointer; }
  .legend-item:hover { text-decoration: underline; }
  .hint { color: var(--muted); font-size: 11px; margin-top: 6px; }
  #wrap { flex: 1; overflow: auto; position: relative; min-height: 200px; }
  svg { display: block; font-family: inherit; }
  .grp, .hit { cursor: pointer; }
  .grp:hover polygon { filter: brightness(0.93); }
  .hit { fill: transparent; }
  .hit:hover { fill: rgba(0, 0, 0, 0.12); }
  #tip {
    position: fixed; pointer-events: none; z-index: 10; display: none;
    background: #fff; border: 1px solid var(--border); border-radius: 6px;
    box-shadow: 0 2px 10px rgba(0,0,0,.12); padding: 6px 9px; font-size: 12px; line-height: 1.45; max-width: 360px;
  }
  #tip b { font-weight: 600; }
  #tip .m { color: var(--muted); }
  @media (max-width: 560px) {
    header { padding: 8px 16px 6px; }
    h1 { font-size: 15px; }
    input[type=search] { width: 190px; }
  }
</style>
</head>
<body>
<div class="banner" id="banner"></div>
<header>
  <h1 id="title"></h1>
  <div class="stats" id="stats"></div>
  <div class="finders">
    <label class="finder"><span>Find sub-clade</span>
      <span class="row"><input type="search" id="find-clade" list="clade-list" placeholder="Type or pick, e.g. AF1a" autocomplete="off">
      <span class="msg" id="clade-msg"></span></span>
    </label>
    <datalist id="clade-list"></datalist>
    <label class="finder"><span>Find accession</span>
      <span class="row"><input type="search" id="find-acc" list="acc-list" placeholder="e.g. AB041966" autocomplete="off">
      <span class="msg" id="acc-msg"></span></span>
    </label>
    <datalist id="acc-list"></datalist>
    <button id="clear" type="button">Clear highlight</button>
  </div>
  <div class="controls">
    <button id="view-sub" type="button">Collapse to sub-clades</button>
    <button id="view-major" type="button">Collapse to major clades</button>
    <button id="expand-all" type="button">Expand all</button>
    <label>Width <input type="range" id="zoom" min="1" max="6" step="0.25" value="1"></label>
    <label>Rows <input type="range" id="row" min="8" max="26" step="1" value="15"></label>
    <button id="download" type="button">Download SVG</button>
  </div>
  <div id="card"></div>
  <details class="legend" open>
    <summary>Major clades (n sequences) &mdash; click to find</summary>
    <div class="legend-items" id="legend"></div>
  </details>
  <div class="hint" id="hint"></div>
</header>
<div id="wrap"></div>
<div id="tip"></div>

<script type="application/json" id="tree-data">__DATA__</script>
<script>
(function () {
  "use strict";

  const D = JSON.parse(document.getElementById("tree-data").textContent);
  const N = D.nodes, SC = D.subclades, MJ = D.majors, META = D.meta, MV = D.metaValues, LIN = D.lin;
  const C = { surface: "#fcfcfb", text: "#0b0b0b", text2: "#52514e", muted: "#8a8984", hl: "#c2185b", band: "#fff1c2", bandStrong: "#ffd666" };
  const FONT = '"DejaVu Sans", Verdana, "Segoe UI", sans-serif';
  const BIG = 4000;   // ask before drawing more sequences than this at once
  const esc = s => String(s).replace(/[&<>"]/g, ch => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[ch]));
  const f1 = v => Math.round(v * 10) / 10;
  const majorName = i => i >= 0 ? MJ[i].name : "mixed";
  const majorColor = i => i >= 0 ? MJ[i].color : D.unknown;
  const minorName = i => { if (i < 0) return "mixed"; const s = SC[i]; return s.major >= 0 && MJ[s.major].name === s.name ? "—" : s.name; };
  const lineageOf = n => n.li ? LIN[n.li] : "not assigned";
  const metaOf = (n, f) => n.m ? MV[f][n.m[f]] : "";
  const metaLabel = c => c.replace(/^collection_/, "").replace(/_/g, " ");

  // ---------------------------------------------------------------
  // Static tree properties (preorder: parent < child)
  // ---------------------------------------------------------------
  const parent = new Int32Array(N.length).fill(-1);
  const X = new Float64Array(N.length), depth = new Float64Array(N.length);
  for (let i = 0; i < N.length; i++) for (const c of N[i].c) { parent[c] = i; X[c] = X[i] + N[c].b; }
  for (let i = N.length - 1; i >= 0; i--) for (const c of N[i].c) depth[i] = Math.max(depth[i], N[c].b + depth[c]);
  let maxX = 0; for (let i = 0; i < N.length; i++) maxX = Math.max(maxX, X[i]);

  function tipsUnder(i) {
    const tips = [], stack = [i];
    while (stack.length) {
      const j = stack.pop();
      if (!N[j].c.length) tips.push(j); else for (let k = N[j].c.length - 1; k >= 0; k--) stack.push(N[j].c[k]);
    }
    return tips;
  }

  const subCollapse = () => { const s = new Set(); N.forEach((n, i) => { if (n.pt) s.add(i); }); return s; };
  const majorCollapse = () => { const s = subCollapse(); N.forEach((n, i) => { if (n.pm) s.add(i); }); return s; };
  let collapsed = subCollapse();

  // Highlight state: highlighted tips below each node
  let hlCount = new Int32Array(N.length), hlTotal = 0, focusTip = -1, pathNodes = new Set();
  function setHighlight(tips) {
    hlCount = new Int32Array(N.length); hlTotal = tips.length;
    for (const t of tips) hlCount[t] = 1;
    for (let i = N.length - 1; i > 0; i--) hlCount[parent[i]] += hlCount[i];
  }
  const isHlTip = i => hlCount[i] > 0 && !N[i].c.length;
  const allHl = i => hlTotal > 0 && hlCount[i] === N[i].nt;
  function reveal(i) { for (let j = parent[i]; j >= 0; j = parent[j]) collapsed.delete(j); }

  // ---------------------------------------------------------------
  // Header, search lists, legend
  // ---------------------------------------------------------------
  document.getElementById("banner").textContent = D.banner;
  document.getElementById("title").textContent = D.title;
  const nTips = N.length ? N[0].nt : 0;
  const realSubs = SC.filter((s, i) => minorName(i) !== "—").length;
  document.getElementById("stats").textContent =
    `${nTips.toLocaleString()} sequences · ${MJ.length} major clades · ${realSubs} minor clades · ` +
    `${D.nLineages.toLocaleString()} with a MADDOG-V2 lineage`;
  document.getElementById("hint").textContent =
    `Click a triangle to expand it (sub-clade → lineage groups → sequences); click a branch point to collapse it. ` +
    `A group is named when at least ${Math.round(D.purity * 100)}% of its sequences share the name; a % shows when it is less than 100%.`;

  const MAJOR_TAG = " (major clade)";
  const cladeOptions = [];
  for (let s = 0; s < SC.length; s++) cladeOptions.push({ value: SC[s].name, label: `${SC[s].major >= 0 ? MJ[SC[s].major].name : ""} · ${SC[s].n} seqs`, sub: s });
  for (let m = 0; m < MJ.length; m++) cladeOptions.push({ value: MJ[m].name + MAJOR_TAG, label: `${MJ[m].n} seqs`, major: m });
  cladeOptions.sort((a, b) => a.value.localeCompare(b.value));
  document.getElementById("clade-list").innerHTML = cladeOptions.map(o => `<option value="${esc(o.value)}" label="${esc(o.label)}"></option>`).join("");

  const accIndex = new Map(), accBase = new Map();
  for (let i = 0; i < N.length; i++) if (!N[i].c.length) {
    const a = N[i].n.toUpperCase();
    accIndex.set(a, i); accBase.set(a.replace(/\.\d+$/, ""), i);
  }
  const accSorted = [...accIndex.keys()].sort();
  const accList = document.getElementById("acc-list");
  function suggestAccessions(prefix) {
    // Only the first matches go in the list: a datalist with every accession is too slow.
    const p = prefix.trim().toUpperCase(), out = [];
    if (p.length >= 2) {
      let lo = 0, hi = accSorted.length;
      while (lo < hi) { const mid = (lo + hi) >> 1; if (accSorted[mid] < p) lo = mid + 1; else hi = mid; }
      for (let j = lo; j < accSorted.length && out.length < 40 && accSorted[j].startsWith(p); j++) out.push(N[accIndex.get(accSorted[j])].n);
    }
    accList.innerHTML = out.map(a => `<option value="${esc(a)}"></option>`).join("");
  }

  const legendEl = document.getElementById("legend");
  legendEl.innerHTML = MJ.map((m, j) =>
    `<span class="legend-item" data-m="${j}"><svg width="22" height="12"><line x1="1" y1="6" x2="21" y2="6" stroke="${m.color}" stroke-width="3"/></svg>${esc(m.name)} (${m.n.toLocaleString()})</span>`).join("");

  // ---------------------------------------------------------------
  // Layout and drawing
  // ---------------------------------------------------------------
  const wrap = document.getElementById("wrap");
  const ui = { zoom: document.getElementById("zoom"), row: document.getElementById("row") };
  const Y = new Float64Array(N.length);
  let geom = null, lastLayout = null;

  function layout() {
    const visible = [], leaves = [];
    let row = 0;
    const stack = [[0, false]];
    while (stack.length) {
      const [i, done] = stack.pop();
      const c = N[i].c;
      if (done) { Y[i] = (Y[c[0]] + Y[c[c.length - 1]]) / 2; continue; }
      visible.push(i);
      if (!c.length || collapsed.has(i)) {
        const h = c.length ? 1.7 : 1.0;
        Y[i] = row + h / 2; row += h; leaves.push(i);
        continue;
      }
      stack.push([i, true]);
      for (let k = c.length - 1; k >= 0; k--) stack.push([c[k], false]);
    }
    return { visible, leaves, rows: row };
  }

  function groupName(n) {
    if (!n.pt) return majorName(n.mj);
    const sub = n.ps !== undefined ? SC[n.ps].name : majorName(n.mj);
    return n.lin ? `${sub} · ${n.lin}` : sub;
  }

  function render() {
    if (!N.length) return;
    const L = layout();
    lastLayout = L;
    const W = wrap.clientWidth, rowH = +ui.row.value, zoom = +ui.zoom.value;
    const colW = [W >= 640 ? 120 : 92, W >= 640 ? 90 : 70, W >= 640 ? 130 : 96];
    const colTotal = colW[0] + colW[1] + colW[2];
    const labelRoom = Math.min(280, Math.max(140, W * 0.22));
    const fs = Math.max(8, Math.min(12, rowH * 0.72));
    const triangleMax = maxX * 0.1;
    const treePx = Math.max(180, W - colTotal - labelRoom - 24) * zoom;
    const k = treePx / (maxX + triangleMax);
    const left = 10, top = 30;
    const svgW = Math.max(W, left + treePx + labelRoom + colTotal + 12);
    const colX = [svgW - colTotal, svgW - colTotal + colW[0], svgW - colTotal + colW[0] + colW[1]];
    const svgH = top + L.rows * rowH + 56;
    const px = x => left + x * k, py = y => top + y * rowH;
    const bands = [], lines = [], hlLines = [], shapes = [], texts = [], hits = [];
    const strong = i => allHl(i) || pathNodes.has(i);
    const branchW = i => strong(i) ? 3 : 1.2;
    const branchC = i => strong(i) ? C.hl : majorColor(N[i].mj);
    const cols = (y, a, b, c, bold, color) => {
      const w = bold ? "700" : "400";
      texts.push(`<text x="${f1(colX[0])}" y="${f1(y)}" font-size="${fs}" font-weight="${w}" fill="${color}" dominant-baseline="central">${esc(a)}</text>`);
      texts.push(`<text x="${f1(colX[1])}" y="${f1(y)}" font-size="${fs}" font-weight="${w}" fill="${bold ? C.text : C.text2}" dominant-baseline="central">${esc(b)}</text>`);
      if (c) texts.push(`<text x="${f1(colX[2])}" y="${f1(y)}" font-size="${fs}" font-weight="${w}" fill="${bold ? C.text : C.muted}" dominant-baseline="central">${esc(c)}</text>`);
    };

    for (const i of L.visible) {
      const n = N[i], x = px(X[i]), y = py(Y[i]);

      if (!n.c.length) {
        const hl = isHlTip(i), focus = i === focusTip;
        if (hl || focus) bands.push(`<rect x="0" y="${f1(y - rowH / 2)}" width="${svgW}" height="${rowH}" fill="${focus ? C.bandStrong : C.band}"/>`);
        shapes.push(`<circle data-i="${i}" cx="${f1(x)}" cy="${f1(y)}" r="${focus ? 5.5 : (hl ? 4 : 2.6)}" fill="${focus ? C.hl : majorColor(n.mj)}" stroke="${C.surface}" stroke-width="1"/>`);
        texts.push(`<text data-i="${i}" x="${f1(x + (focus ? 9 : 6))}" y="${f1(y)}" font-size="${fs}" font-weight="${hl || focus ? 700 : 400}" fill="${hl || focus ? C.text : C.text2}" dominant-baseline="central">${esc(n.n)}</text>`);
        cols(y, majorName(n.mj), minorName(n.sc), n.li ? LIN[n.li] : "", hl || focus, majorColor(n.mj));
        continue;
      }

      if (collapsed.has(i)) {
        const h = hlCount[i], all = allHl(i);
        const d = Math.min(Math.max(depth[i], maxX * 0.012), triangleMax) * k, half = 0.62 * rowH;
        if (h) bands.push(`<rect x="0" y="${f1(y - 0.85 * rowH)}" width="${svgW}" height="${f1(1.7 * rowH)}" fill="${C.band}"/>`);
        const label = groupName(n) + (n.dp ? ` (${n.dp}%)` : "") + `  ·  ${n.nt.toLocaleString()} seqs` +
          (h && !all ? `  ·  ${h.toLocaleString()} highlighted` : "");
        shapes.push(`<g class="grp" data-i="${i}" data-act="expand"><polygon points="${f1(x)},${f1(y)} ${f1(x + d)},${f1(y - half)} ${f1(x + d)},${f1(y + half)}" fill="${majorColor(n.mj)}" fill-opacity="${h ? 0.85 : 0.35}" stroke="${h ? C.hl : majorColor(n.mj)}" stroke-width="${h ? 2 : 0.8}"/>` +
          `<text x="${f1(x + d + 6)}" y="${f1(y)}" font-size="${fs + 1}" font-weight="700" fill="${C.text}" dominant-baseline="central">${esc(label)}</text></g>`);
        const minor = n.ps !== undefined ? minorName(n.ps) : minorName(n.sc);
        cols(y, majorName(n.mj), minor, n.lin || "", true, majorColor(n.mj));
        continue;
      }

      const ys = n.c.map(c => py(Y[c]));
      (strong(i) ? hlLines : lines).push(`<line x1="${f1(x)}" y1="${f1(Math.min(...ys))}" x2="${f1(x)}" y2="${f1(Math.max(...ys))}" stroke="${branchC(i)}" stroke-width="${branchW(i)}"/>`);
      for (const c of n.c) {
        const yc = py(Y[c]);
        (strong(c) ? hlLines : lines).push(`<line x1="${f1(x)}" y1="${f1(yc)}" x2="${f1(px(X[c]))}" y2="${f1(yc)}" stroke="${branchC(c)}" stroke-width="${branchW(c)}"/>`);
      }
      if (i !== 0) hits.push(`<circle class="hit" data-i="${i}" data-act="collapse" cx="${f1(x)}" cy="${f1(y)}" r="6"/>`);
    }

    const length = Math.pow(10, Math.floor(Math.log10(maxX * 0.2)));
    const sbX = left + 4, sbY = svgH - 22;
    const scale = `<line x1="${f1(sbX)}" y1="${sbY}" x2="${f1(sbX + length * k)}" y2="${sbY}" stroke="${C.text2}" stroke-width="1.2"/>` +
      `<text x="${f1(sbX)}" y="${sbY - 5}" font-size="11" fill="${C.text2}">${+length.toPrecision(3)} subs/site</text>`;
    const head = (x, label) => `<text x="${f1(x)}" y="14" font-size="12" font-weight="700" fill="${C.text}" dominant-baseline="central">${esc(label)}</text>`;

    wrap.innerHTML =
      `<svg id="tree" xmlns="http://www.w3.org/2000/svg" width="${f1(svgW)}" height="${f1(svgH)}" font-family='${FONT}'>` +
      `<rect width="100%" height="100%" fill="${C.surface}"/>` + bands.join("") +
      `<g stroke-linecap="round">${lines.join("")}${hlLines.join("")}${shapes.join("")}${texts.join("")}${hits.join("")}</g>` +
      head(colX[0], "Major clade") + head(colX[1], "Minor clade") + head(colX[2], "Lineage") + scale + `</svg>`;
    geom = { py };
  }

  // Sequences that would be drawn if node i were expanded one step (for the size warning).
  function shownTips(set) {
    let count = 0;
    const stack = [0];
    while (stack.length) {
      const i = stack.pop();
      if (!N[i].c.length) count++;
      else if (!set.has(i)) for (const c of N[i].c) stack.push(c);
    }
    return count;
  }
  function okToDraw(set) {
    const n = shownTips(set);
    return n <= BIG || confirm(`This draws ${n.toLocaleString()} sequences and may be slow. Continue?`);
  }

  // ---------------------------------------------------------------
  // Interaction
  // ---------------------------------------------------------------
  function scrollToNode(i) {
    if (geom) wrap.scrollTo({ top: Math.max(0, geom.py(Y[i]) - wrap.clientHeight / 3), behavior: "smooth" });
  }
  // Scroll to the drawn row that holds the most highlighted sequences.
  function scrollToHighlight() {
    if (!lastLayout || !hlTotal) return;
    let best = -1;
    for (const i of lastLayout.leaves) if (hlCount[i] && (best < 0 || hlCount[i] > hlCount[best])) best = i;
    if (best >= 0) scrollToNode(best);
  }
  const card = document.getElementById("card");
  const cladeMsg = document.getElementById("clade-msg"), accMsg = document.getElementById("acc-msg");
  const findClade = document.getElementById("find-clade"), findAcc = document.getElementById("find-acc");

  function clearAll() {
    setHighlight([]); focusTip = -1; pathNodes = new Set();
    card.style.display = "none"; cladeMsg.textContent = ""; accMsg.textContent = "";
  }

  function selectClade(value) {
    const v = value.trim();
    if (!v) return;
    const opt = cladeOptions.find(o => o.value.toLowerCase() === v.toLowerCase()) ||
                cladeOptions.find(o => o.value.toLowerCase() === (v + MAJOR_TAG).toLowerCase());
    if (!opt) { cladeMsg.className = "msg err"; cladeMsg.textContent = "no such sub-clade"; return; }
    focusTip = -1; pathNodes = new Set(); card.style.display = "none"; accMsg.textContent = "";
    const tips = [];
    for (let i = 0; i < N.length; i++) if (!N[i].c.length &&
        (opt.major !== undefined ? N[i].mj === opt.major : N[i].sc === opt.sub)) tips.push(i);
    // Back to the sub-clade view: the clade's triangles light up, and triangles of other
    // clades that hold some of its sequences show how many.
    collapsed = subCollapse();
    setHighlight(tips);
    let groups = 0;
    render();
    for (const i of lastLayout.leaves) if (hlCount[i]) groups++;
    cladeMsg.className = "msg";
    cladeMsg.textContent = `${tips.length.toLocaleString()} sequence${tips.length === 1 ? "" : "s"} highlighted in ${groups} row${groups === 1 ? "" : "s"}`;
    scrollToHighlight();
  }

  function showCard(i) {
    const n = N[i];
    const fields = [["Accession", n.n], ["Major clade", majorName(n.mj)], ["Minor clade", minorName(n.sc)], ["Lineage", lineageOf(n)]];
    if (n.lp) fields.push(["Parent lineage", LIN[n.lp]]);
    META.forEach((c, f) => { const v = metaOf(n, f); if (v) fields.push([metaLabel(c), v]); });
    card.innerHTML = `<button class="close" type="button" aria-label="Close">×</button><div class="grid">` +
      fields.map(([k, v], j) => `<div><div class="k">${esc(k)}</div><div class="v"${j === 1 ? ` style="color:${majorColor(n.mj)}"` : ""}>${esc(v)}</div></div>`).join("") + `</div>`;
    card.style.display = "block";
  }

  function selectAccession(value) {
    const v = value.trim().toUpperCase();
    if (!v) return;
    const i = accIndex.has(v) ? accIndex.get(v) : accBase.get(v.replace(/\.\d+$/, ""));
    if (i === undefined) {
      accMsg.className = "msg err"; accMsg.textContent = "not in the tree"; card.style.display = "none"; return;
    }
    accMsg.className = "msg"; accMsg.textContent = ""; cladeMsg.textContent = "";
    focusTip = i; pathNodes = new Set();
    for (let j = i; j >= 0; j = parent[j]) pathNodes.add(j);
    showCard(i);
    collapsed = subCollapse();
    reveal(i);
    setHighlight([i]);
    render(); scrollToNode(i);
  }

  findClade.addEventListener("change", () => selectClade(findClade.value));
  findClade.addEventListener("keydown", ev => { if (ev.key === "Enter") { ev.preventDefault(); selectClade(findClade.value); } });
  findClade.addEventListener("input", () => {
    if (cladeOptions.some(o => o.value === findClade.value)) selectClade(findClade.value);
    else if (!findClade.value) cladeMsg.textContent = "";
  });
  findAcc.addEventListener("change", () => selectAccession(findAcc.value));
  findAcc.addEventListener("keydown", ev => { if (ev.key === "Enter") { ev.preventDefault(); selectAccession(findAcc.value); } });
  findAcc.addEventListener("input", () => {
    suggestAccessions(findAcc.value);
    if (accIndex.has(findAcc.value.trim().toUpperCase()) && accList.options.length === 1) selectAccession(findAcc.value);
    else if (!findAcc.value) accMsg.textContent = "";
  });
  card.addEventListener("click", ev => { if (ev.target.closest(".close")) { clearAll(); findAcc.value = ""; render(); } });

  legendEl.addEventListener("click", ev => {
    const el = ev.target.closest("[data-m]");
    if (!el) return;
    findClade.value = MJ[+el.dataset.m].name + MAJOR_TAG;
    selectClade(findClade.value);
  });

  wrap.addEventListener("click", ev => {
    const el = ev.target.closest("[data-act]");
    if (!el) return;
    const i = +el.dataset.i;
    if (el.dataset.act === "expand") {
      const next = new Set(collapsed); next.delete(i);
      if (!okToDraw(next)) return;
      collapsed = next;
    } else collapsed.add(i);
    render();
  });

  const tip = document.getElementById("tip");
  wrap.addEventListener("mousemove", ev => {
    const el = ev.target.closest("[data-i]");
    if (!el) { tip.style.display = "none"; return; }
    const i = +el.dataset.i, n = N[i];
    let body;
    if (!n.c.length) {
      body = `<b>${esc(n.n)}</b><br><span class="m">major clade</span> <b style="color:${majorColor(n.mj)}">${esc(majorName(n.mj))}</b>` +
        `<br><span class="m">minor clade</span> ${esc(minorName(n.sc))}` +
        `<br><span class="m">lineage</span> ${esc(lineageOf(n))}` + (n.lp ? ` <span class="m">(parent ${esc(LIN[n.lp])})</span>` : "") +
        META.map((c, f) => { const v = metaOf(n, f); return v ? `<br><span class="m">${esc(metaLabel(c))}</span> ${esc(v)}` : ""; }).join("");
    } else {
      const counts = new Map();
      for (const t of tipsUnder(i)) { const s = N[t].sc; counts.set(s, (counts.get(s) || 0) + 1); }
      const ranked = [...counts].sort((a, b) => b[1] - a[1]);
      body = `<b>${esc(n.pt || n.pm ? groupName(n) : n.nt.toLocaleString() + " sequences")}</b>` +
        (n.pt || n.pm ? `<br>${n.nt.toLocaleString()} sequences` + (n.dp ? ` · ${n.dp}% share the name` : "") : "") +
        `<br><span class="m">sub-clades</span> ${ranked.slice(0, 6).map(([s, c]) => `${esc(s >= 0 ? SC[s].name : "?")} (${c})`).join(", ")}${ranked.length > 6 ? ` +${ranked.length - 6}` : ""}` +
        (n.s ? `<br><span class="m">support</span> ${esc(n.s)}` : "") +
        `<br><span class="m">${collapsed.has(i) ? "Click to expand" : "Click to collapse"}</span>`;
    }
    tip.innerHTML = body;
    tip.style.display = "block";
    tip.style.left = Math.min(ev.clientX + 14, window.innerWidth - tip.offsetWidth - 8) + "px";
    tip.style.top = Math.min(ev.clientY + 14, window.innerHeight - tip.offsetHeight - 8) + "px";
  });
  wrap.addEventListener("mouseleave", () => { tip.style.display = "none"; });

  document.getElementById("clear").addEventListener("click", () => { clearAll(); findClade.value = ""; findAcc.value = ""; render(); });
  document.getElementById("view-sub").addEventListener("click", () => { collapsed = subCollapse(); if (focusTip >= 0) reveal(focusTip); render(); scrollToHighlight(); });
  document.getElementById("view-major").addEventListener("click", () => { collapsed = majorCollapse(); render(); scrollToHighlight(); });
  document.getElementById("expand-all").addEventListener("click", () => { if (okToDraw(new Set())) { collapsed = new Set(); render(); } });
  for (const el of [ui.zoom, ui.row]) el.addEventListener("input", render);

  document.getElementById("download").addEventListener("click", () => {
    const svg = document.getElementById("tree");
    if (!svg) return;
    const blob = new Blob([new XMLSerializer().serializeToString(svg)], { type: "image/svg+xml" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = "rabv_tree.svg"; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  });

  let lastWidth = 0;
  new ResizeObserver(() => { if (wrap.clientWidth !== lastWidth) { lastWidth = wrap.clientWidth; render(); } }).observe(wrap);
})();
</script>
</body>
</html>
"""


def parse_args():
    ref = HERE / "reference" / "rabv"
    default_db = HERE / "rabv-gDB_16062026.db"
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree", default=str(HERE / "extract_sub_tree" / "full_tree_midpoint_rooted.contree"),
                        help="Newick tree to draw.")
    parser.add_argument("--references", default=None,
                        help="Draw only these tips (a Newick tree or one accession per line), "
                             "e.g. reference/rabv/tree/ref.treefile. Default: every tip of --tree.")
    parser.add_argument("--db", default=str(default_db) if default_db.exists() else None,
                        help="rabv-gDB SQLite file: EPA major/minor clade, country, host, year.")
    parser.add_argument("--major", default=str(ref / "epa_clades" / "ref_major_clades.tsv"),
                        help="TSV sequence_id, major_clade (used before the database).")
    parser.add_argument("--minor", default=str(ref / "epa_clades" / "ref_minor_clades.tsv"),
                        help="TSV sequence_id, minor_clade (NULL = none).")
    parser.add_argument("--lineage_dir", default=str(HERE / "extract_sub_tree" / "subclades"),
                        help="Directory with <subclade>/meta_data.tsv (ID, Lineage, Parent) or "
                             "<subclade>/<subclade>.tip_lineage_assignments.tsv.")
    parser.add_argument("--metadata", default=None,
                        help="Optional TSV/CSV with an accession column and country/host/collection_year.")
    parser.add_argument("--purity", type=float, default=0.9,
                        help="Share of a branch's sequences that must have one name to collapse it "
                             "under that name (default 0.9).")
    parser.add_argument("--outgroups", nargs="+", default=OUTGROUPS,
                        help="Root on these tips when they are in the tree.")
    parser.add_argument("--keep_outgroups", action="store_true", help="Draw the outgroups instead of dropping them.")
    parser.add_argument("--title", default="RABV phylogeny: sub-clades and MADDOG-V2 lineages", help="Page title.")
    parser.add_argument("--output", default="RABV_reference_tree.html", help="Output .html path.")
    return parser.parse_args()
TREE_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/rabv_tree.nwk"
METADATA_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/rabv_meta_data.csv"
OUTPUT_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/treeData.json"
DB = "/Volumes/My Passport/CVR/gdb/RABV/rabv-gDB_11112026.db"
def main():
    # args = parse_args()
    # tree = Phylo.read(args.tree, "newick")
    tree = Phylo.read(TREE_FILE, "newick")
    # if args.references:
    #     keep = tip_names(None)
    #     before = tree.count_terminals()
    #     tree = subset_tree(tree, keep)
    #     print(f"Reduced {args.tree} from {before} to {tree.count_terminals()} tips in {args.references}")
    # tree = root_tree(tree, OUTGROUPS, False)
    tips = TipIndex([t.name for t in tree.get_terminals() if t.name])

    clades = read_clades(tips, None, None, DB)
    tip_major, tip_sub = {}, {}
    for tip in tips.exact:
        major, minor = clades.get(tip, ("", ""))
        major = major or UNCLASSIFIED
        tip_major[tip] = major
        tip_sub[tip] = minor if minor and major != UNCLASSIFIED else major
    unclassified = sum(1 for v in tip_major.values() if v == UNCLASSIFIED)
    if unclassified:
        print(f"Note: {unclassified} tips have no EPA major clade; shown as '{UNCLASSIFIED}'")
    lineage_dir = False
    lineages = read_lineages(lineage_dir, tip_sub, tips) if lineage_dir else {}

    meta_values = read_db_metadata(DB, tips) if DB else {c: {} for c in META_COLUMNS}
    if METADATA_FILE:
        found, values = read_metadata(METADATA_FILE, META_COLUMNS, None, tips)
        for column in found:
            for tip, value in values[column].items():
                meta_values[column].setdefault(tip, value)

    nodes, end = flat_nodes(tree)
    subclades, majors = annotate(nodes, end, tip_major, tip_sub, lineages, 0.9)
    lin_list, meta_lists = encode_tips(nodes, lineages, meta_values)
    n_lineages = sum(1 for v in lineages.values() if v[0] and v[0] != "Unassigned")
    data = {
        "title": "test", "banner": BANNER, "purity": 0.9, "unknown": UNKNOWN,
        "nLineages": n_lineages, "meta": META_COLUMNS, "metaValues": meta_lists, "lin": lin_list,
        "subclades": subclades, "majors": majors, "nodes": nodes,
    }
    # path = write_html(args.output, data)
    # print(f"Tips: {len(tips.exact)}  major clades: {len(majors)}  sub-clades: {len(subclades)}  "
    #       f"with lineage: {n_lineages}")
    # print(f"Wrote {path}  ({path.stat().st_size / 1e6:.1f} MB)")

    with open(OUTPUT_FILE, "w", encoding="utf-8") as file:
    
        json.dump(
            data,
            file,
            indent=2
        )
    return 0


if __name__ == "__main__":
    main()
    # sys.setrecursionlimit(100_000)
    # raise SystemExit(main())
