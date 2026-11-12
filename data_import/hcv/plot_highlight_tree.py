from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
import threading
from pathlib import Path
import json
from Bio import Phylo


PALETTE = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
    "#7cb6f0", "#a33a12", "#5fc9a0", "#9c6b00",
    "#b0437a", "#7fbf3f", "#9d8fe0", "#8f1f1f",
]
OTHER = "#6f6e69"          # groups beyond the palette
BACKGROUND = "#b9b8b3"     # branches without highlighted tips
MIXED = "#52514e"          # branches above highlighted tips of several groups
GROUP_FILL = "#e4e3df"
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#8a8984"

ID_COLUMNS = ["accession", "primary_accession", "ID", "sequence_id", "accession_version"]
HEADER_WORDS = {c.lower() for c in ID_COLUMNS} | {"id", "acc"}
VERSION = re.compile(r"\.\d+$")


# ------------------------------------------------------------
# Tips: accessions match with or without a version (GU798494 = GU798494.1)
# ------------------------------------------------------------

class TipIndex:
    def __init__(self, names):
        self.exact = set(names)
        self.base = {}
        for name in names:
            self.base.setdefault(VERSION.sub("", name), name)

    def find(self, accession):
        accession = accession.strip()
        if accession in self.exact:
            return accession
        return self.base.get(VERSION.sub("", accession))


# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

def open_table(path):
    path = Path(path)
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    csv.field_size_limit(sys.maxsize)
    handle = path.open(newline="", encoding="utf-8")
    return handle, csv.DictReader(handle, delimiter=delimiter)


def read_metadata(path, columns, id_column, tips):
    """Return (found columns, {column: {tip: value}})."""
    handle, reader = open_table(path)
    with handle:
        fields = reader.fieldnames or []
        by_lower = {f.lower(): f for f in fields}

        if id_column:
            if id_column not in fields:
                raise ValueError(f"--id_column {id_column} not in {path}")
        else:
            id_column = next((by_lower[c.lower()] for c in ID_COLUMNS if c.lower() in by_lower), None)
            if id_column is None:
                raise ValueError(
                    f"{path} has no accession column ({', '.join(ID_COLUMNS)}); "
                    "set --id_column"
                )

        # Exact name first, else a prefixed one (year -> collection_year).
        sources = {}
        for column in columns:
            source = by_lower.get(column.lower()) or next(
                (f for f in fields if f.lower().endswith("_" + column.lower())), None
            )
            if source is None:
                print(f"Warning: column '{column}' not found in {path}")
            else:
                sources[column] = source

        values = {column: {} for column in sources}
        matched = 0
        for row in reader:
            tip = tips.find(row.get(id_column) or "")
            if tip is None:
                continue
            matched += 1
            for column, source in sources.items():
                value = (row.get(source) or "").strip()
                if value and tip not in values[column]:
                    values[column][tip] = value

    print(f"Metadata: {matched} rows matched tree tips; columns {list(sources)}")
    return list(sources), values


def read_highlights(path, tips):
    """Return ({tip: group}, [groups in file order])."""
    tip_group, groups, missing = {}, [], []
    with open(path, encoding="utf-8") as handle:
        for n, line in enumerate(handle):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in re.split(r"\t|,", line) if p.strip()]
            if len(parts) == 1 and len(line.split()) > 1:
                parts = line.split()
            accession = parts[0]
            group = parts[1] if len(parts) > 1 else "Highlighted"
            tip = tips.find(accession)
            if tip is None:
                if n == 0 and accession.lower() in HEADER_WORDS:
                    continue   # header line
                missing.append(accession)
                continue
            if tip not in tip_group:
                tip_group[tip] = group
                if group not in groups:
                    groups.append(group)

    print(f"Highlight: {len(tip_group)} accessions in the tree, {len(groups)} groups")
    if missing:
        shown = ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
        print(f"Warning: {len(missing)} highlight accessions not in the tree: {shown}")
    return tip_group, groups


# ------------------------------------------------------------
# Tree -> flat preorder node list for the page (iterative: deep trees are fine)
# ------------------------------------------------------------

def tree_nodes(tree, tip_group, groups, meta_columns, meta_values):
    group_index = {g: i for i, g in enumerate(groups)}

    order, stack = [], [tree.root]
    while stack:
        clade = stack.pop()
        order.append(clade)
        stack.extend(clade.clades)

    # Tip counts, then ladderize (smaller subtrees first, like Phylo.ladderize).
    count = {}
    for clade in reversed(order):
        count[clade] = 1 if clade.is_terminal() else sum(count[c] for c in clade.clades)
    for clade in order:
        clade.clades.sort(key=lambda c: count[c])

    nodes, index = [], {}
    stack = [(tree.root, -1)]
    while stack:
        clade, parent = stack.pop()
        i = len(nodes)
        index[clade] = i
        node = {"c": [], "b": clade.branch_length or 0.0}
        if clade.is_terminal():
            node["n"] = clade.name or ""
            if node["n"] in tip_group:
                node["g"] = group_index[tip_group[node["n"]]]
            if meta_columns:
                node["m"] = [meta_values[c].get(node["n"], "") for c in meta_columns]
        else:
            support = clade.confidence if clade.confidence is not None else clade.name
            node["s"] = "" if support is None else str(support)
        nodes.append(node)
        if parent >= 0:
            nodes[parent]["c"].append(i)
        for child in reversed(clade.clades):
            stack.append((child, i))

    # Per node: tips, highlighted tips, and the single group below it (-1 none, -2 mixed).
    for node in reversed(nodes):
        if not node["c"]:
            node["nt"] = 1
            node["nh"] = 1 if "g" in node else 0
            node["pg"] = node.get("g", -1) if "g" in node else -1
            node["pure"] = "g" in node
        else:
            kids = [nodes[c] for c in node["c"]]
            node["nt"] = sum(k["nt"] for k in kids)
            node["nh"] = sum(k["nh"] for k in kids)
            node["pure"] = all(k["pure"] for k in kids) and len({k["pg"] for k in kids}) == 1
            node["pg"] = kids[0]["pg"] if node["pure"] else -2
    for node in nodes:
        # Branch colour index: group when every tip below is that group.
        node["k"] = node["pg"] if node["pure"] else (-2 if node["nh"] else -1)
        del node["pure"], node["pg"]

    # Label each group at its MRCA when the group is a clade on its own.
    size = {}
    for g in tip_group.values():
        size[group_index[g]] = size.get(group_index[g], 0) + 1
    labelled = set()
    for node in nodes:   # preorder: the first pure node of full size is the MRCA
        g = node["k"]
        if g >= 0 and node["c"] and g not in labelled and node["nt"] == size[g]:
            node["l"] = g
            labelled.add(g)

    mrca = index[tree.common_ancestor(list(tip_group))] if tip_group else -1
    return nodes, mrca


def group_colors(groups):
    return [PALETTE[i] if i < len(PALETTE) else OTHER for i in range(len(groups))]


def write_html(output, title, nodes, mrca, groups, meta_columns):
    data = {
        "title": title,
        "mrca": mrca,
        "groups": groups,
        "groupColors": group_colors(groups),
        "meta": meta_columns,
        "palette": PALETTE,
        "other": OTHER,
        "nodes": nodes,
    }
    # payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    # page = HTML_TEMPLATE.replace("__TITLE__", html.escape(title)).replace("__DATA__", payload)
    # for key, value in {
    #     "__SURFACE__": SURFACE, "__TEXT_PRIMARY__": TEXT_PRIMARY,
    #     "__TEXT_SECONDARY__": TEXT_SECONDARY, "__TEXT_MUTED__": TEXT_MUTED,
    #     "__BACKGROUND__": BACKGROUND, "__MIXED__": MIXED, "__GROUP_FILL__": GROUP_FILL,
    # }.items():
    #     page = page.replace(key, value)

    # output = Path(output)
    # if output.suffix.lower() != ".html":
    #     output = output.with_suffix(".html")
    # output.parent.mkdir(parents=True, exist_ok=True)
    # output.write_text(page, encoding="utf-8")
    return data


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  :root {
    --surface: __SURFACE__; --text: __TEXT_PRIMARY__; --text-2: __TEXT_SECONDARY__;
    --muted: __TEXT_MUTED__; --border: #e4e3df; --match: #ffe9a8;
  }
  * { box-sizing: border-box; }
  html, body { margin: 0; height: 100%; }
  body {
    background: var(--surface); color: var(--text);
    font-family: "DejaVu Sans", Verdana, "Segoe UI", sans-serif; font-size: 13px;
    display: flex; flex-direction: column;
  }
  header { flex: none; padding: 12px 16px 8px; border-bottom: 1px solid var(--border); }
  h1 { font-size: 17px; font-weight: 600; margin: 0 0 2px; }
  .stats { color: var(--text-2); font-size: 12px; margin-bottom: 8px; }
  .controls { display: flex; flex-wrap: wrap; gap: 6px 14px; align-items: center; }
  #meta-cols { display: inline-flex; flex-wrap: wrap; gap: 4px 10px; }
  .controls label { display: inline-flex; align-items: center; gap: 5px; color: var(--text-2); font-size: 12px; }
  .search { display: inline-flex; align-items: center; gap: 6px; }
  input[type=search] {
    font: inherit; padding: 4px 8px; width: 210px; max-width: 60vw;
    border: 1px solid #cfcec9; border-radius: 5px; background: #fff; color: var(--text);
  }
  #match-count { color: var(--muted); font-size: 12px; min-width: 52px; }
  button, select {
    font: inherit; font-size: 12px; padding: 4px 9px; cursor: pointer;
    border: 1px solid #cfcec9; border-radius: 5px; background: #fff; color: var(--text);
  }
  select { padding: 3px 4px; }
  button:hover { background: #f1f0ed; }
  input[type=range] { width: 90px; }
  details.legend { margin-top: 8px; }
  details.legend summary { cursor: pointer; font-size: 12px; font-weight: 600; }
  .legend-items { display: flex; flex-wrap: wrap; gap: 4px 16px; margin-top: 6px; }
  .legend-item { display: inline-flex; align-items: center; gap: 6px; font-size: 12px; cursor: default; }
  .legend-item.jump { cursor: pointer; }
  .legend-item.jump:hover { text-decoration: underline; }
  .hint { color: var(--muted); font-size: 11px; margin-top: 6px; }
  #wrap { flex: 1; overflow: auto; position: relative; }
  svg { display: block; font-family: inherit; }
  .hit, .grp { cursor: pointer; }
  .grp:hover polygon { fill: #d6d5d0; }
  .hit { fill: transparent; }
  .hit:hover { fill: rgba(0, 0, 0, 0.12); }
  #tip {
    position: fixed; pointer-events: none; z-index: 10; display: none;
    background: #fff; border: 1px solid var(--border); border-radius: 6px;
    box-shadow: 0 2px 10px rgba(0,0,0,.12); padding: 6px 9px; font-size: 12px; line-height: 1.45;
    max-width: 320px;
  }
  #tip b { font-weight: 600; }
  #tip .m { color: var(--muted); }
  @media (max-width: 560px) {
    header { padding: 10px 16px 6px; }
    h1 { font-size: 15px; }
    input[type=search] { width: 150px; }
  }
</style>
</head>
<body>
<header>
  <h1 id="title"></h1>
  <div class="stats" id="stats"></div>
  <div class="controls">
    <span class="search">
      <input type="search" id="search" placeholder="Search accession, group or metadata" autocomplete="off">
      <span id="match-count"></span>
    </span>
    <button id="expand-all" type="button">Expand all</button>
    <button id="reset" type="button">Collapse unhighlighted</button>
    <label><input type="checkbox" id="truncate" checked> Truncate root branch</label>
    <label>Colour tips by <select id="color-by"><option value="-1">Highlight group</option></select></label>
    <span id="meta-cols"></span>
    <label>Width <input type="range" id="zoom" min="1" max="8" step="0.25" value="1"></label>
    <label>Rows <input type="range" id="row" min="8" max="26" step="1" value="14"></label>
    <button id="download" type="button">Download SVG</button>
  </div>
  <details class="legend" open>
    <summary id="legend-title"></summary>
    <div class="legend-items" id="legend"></div>
  </details>
  <div class="hint">Click a grey triangle to expand it; click a branch point to collapse it. Click a group in the legend to jump to it. Enter in the search box jumps to the next match.</div>
</header>
<div id="wrap"></div>
<div id="tip"></div>

<script type="application/json" id="tree-data">__DATA__</script>
<script>
(function () {
  "use strict";

  const D = JSON.parse(document.getElementById("tree-data").textContent);
  const N = D.nodes, G = D.groups, GC = D.groupColors;
  const C = {
    surface: "__SURFACE__", text: "__TEXT_PRIMARY__", text2: "__TEXT_SECONDARY__",
    muted: "__TEXT_MUTED__", background: "__BACKGROUND__", mixed: "__MIXED__", group: "__GROUP_FILL__",
  };
  const FONT = '"DejaVu Sans", Verdana, "Segoe UI", sans-serif';
  const esc = s => String(s).replace(/[&<>"]/g, ch => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[ch]));
  const f1 = v => Math.round(v * 10) / 10;
  const cap = s => s ? (s[0].toUpperCase() + s.slice(1)).replace(/_/g, " ") : s;
  // Shorten a column value to roughly fit its width (full value is in the tooltip).
  const fit = (v, w) => { const n = Math.floor((w - 10) / 6.5); return v.length > n ? v.slice(0, n - 1) + "…" : v; };
  const branchColor =n => n.k >= 0 ? GC[n.k] : (n.k === -2 ? C.mixed : C.background);

  // ---------------------------------------------------------------
  // Static tree properties (nodes are in preorder: parent < child)
  // ---------------------------------------------------------------
  const parent = new Int32Array(N.length).fill(-1);
  const X = new Float64Array(N.length), depth = new Float64Array(N.length);
  for (let i = 0; i < N.length; i++) for (const c of N[i].c) { parent[c] = i; X[c] = X[i] + N[c].b; }
  for (let i = N.length - 1; i >= 0; i--) for (const c of N[i].c) depth[i] = Math.max(depth[i], N[c].b + depth[c]);

  // Subtrees without highlighted tips that hang off highlighted branches start collapsed.
  function defaultCollapsed() {
    const set = new Set();
    if (!N.length || !N[0].nh) return set;
    for (let i = 1; i < N.length; i++) if (N[i].c.length && !N[i].nh && N[parent[i]].nh) set.add(i);
    return set;
  }
  let collapsed = defaultCollapsed();

  const tipCache = new Map();
  function tipsUnder(i) {
    if (tipCache.has(i)) return tipCache.get(i);
    const tips = [], stack = [i];
    while (stack.length) {
      const j = stack.pop();
      if (!N[j].c.length) tips.push(j);
      else for (let k = N[j].c.length - 1; k >= 0; k--) stack.push(N[j].c[k]);
    }
    tipCache.set(i, tips);
    return tips;
  }

  // ---------------------------------------------------------------
  // Metadata colour schemes: numeric (year) -> blue ramp, else palette by frequency
  // ---------------------------------------------------------------
  const META = D.meta || [];
  const metaValue = (i, f) => (N[i].m && N[i].m[f]) || "";
  const schemes = META.map((name, f) => {
    const counts = new Map();
    for (const n of N) { const v = n.m && n.m[f]; if (!n.c.length && v) counts.set(v, (counts.get(v) || 0) + 1); }
    const values = [...counts.keys()];
    if (values.length > 8 && values.every(v => /^-?\d+(\.\d+)?$/.test(v))) {
      const nums = values.map(Number), lo = Math.min(...nums), hi = Math.max(...nums);
      const a = [214, 230, 247], b = [8, 48, 107];
      const color = v => { const t = hi > lo ? (Number(v) - lo) / (hi - lo) : 1;
        return `rgb(${a.map((c, j) => Math.round(c + (b[j] - c) * t)).join(",")})`; };
      const legend = Array.from({ length: 5 }, (_, s) => {
        const v = Math.round(lo + (hi - lo) * s / 4);
        return { label: String(v), color: color(v), kind: "dot" };
      });
      return { color, legend };
    }
    const ranked = values.sort((x, y) => counts.get(y) - counts.get(x) || x.localeCompare(y));
    const slot = new Map(ranked.slice(0, D.palette.length).map((v, j) => [v, D.palette[j]]));
    const legend = ranked.slice(0, D.palette.length).map(v => ({ label: `${v} (${counts.get(v)})`, color: slot.get(v), kind: "dot" }));
    if (ranked.length > D.palette.length) legend.push({ label: `Other (${ranked.length - D.palette.length})`, color: D.other, kind: "dot" });
    return { color: v => slot.get(v) || D.other, legend };
  });

  function topValues(i, f, limit) {
    const counts = new Map();
    for (const t of tipsUnder(i)) { const v = metaValue(t, f); if (v) counts.set(v, (counts.get(v) || 0) + 1); }
    const ranked = [...counts.entries()].sort((a, b) => b[1] - a[1]);
    return ranked.slice(0, limit).map(([v, c]) => `${esc(v)} (${c})`).join(", ") +
      (ranked.length > limit ? ` +${ranked.length - limit}` : "");
  }

  function collapsedLabel(i) {
    const top = META.length ? topValues(i, 0, 2) : "";
    return `${N[i].nt} seqs` + (top ? `  ·  ${top.replace(/&amp;/g, "&")}` : "");
  }

  // ---------------------------------------------------------------
  // Header, controls, legend
  // ---------------------------------------------------------------
  document.getElementById("title").textContent = D.title;
  const groupSize = G.map(() => 0);
  for (const n of N) if (!n.c.length && n.g !== undefined) groupSize[n.g]++;
  document.getElementById("stats").textContent =
    `${N.length ? N[0].nt : 0} sequences · ${N.length ? N[0].nh : 0} highlighted · ${G.length} groups`;

  const wrap = document.getElementById("wrap");
  const ui = {
    search: document.getElementById("search"), count: document.getElementById("match-count"),
    truncate: document.getElementById("truncate"), zoom: document.getElementById("zoom"),
    row: document.getElementById("row"), colorBy: document.getElementById("color-by"),
    metaCols: document.getElementById("meta-cols"),
  };
  if (D.mrca < 0) { ui.truncate.checked = false; ui.truncate.disabled = true; }
  ui.colorBy.innerHTML += META.map((name, f) => `<option value="${f}">${esc(name)}</option>`).join("");
  ui.metaCols.innerHTML = META.map((name, f) =>
    `<label><input type="checkbox" data-f="${f}" checked> ${esc(cap(name))}</label>`).join("");

  function mark(item) {
    if (item.kind === "line") return `<svg width="22" height="12"><line x1="1" y1="6" x2="21" y2="6" stroke="${item.color}" stroke-width="2"/><circle cx="11" cy="6" r="3.5" fill="${item.color}"/></svg>`;
    if (item.kind === "plain") return `<svg width="22" height="12"><line x1="1" y1="6" x2="21" y2="6" stroke="${item.color}" stroke-width="2"/></svg>`;
    if (item.kind === "patch") return `<svg width="18" height="12"><rect x="1" y="1" width="16" height="10" fill="${item.color}" stroke="${C.background}" stroke-width="0.8"/></svg>`;
    return `<svg width="12" height="12"><circle cx="6" cy="6" r="4.5" fill="${item.color}"/></svg>`;
  }
  const legendEl = document.getElementById("legend"), legendTitle = document.getElementById("legend-title");
  function updateLegend() {
    const f = +ui.colorBy.value;
    let items;
    if (f < 0) {
      legendTitle.textContent = "Highlight group (n sequences)";
      items = G.map((g, j) => ({ label: `${g} (${groupSize[j]})`, color: GC[j], kind: "line", g: j }));
      items.push({ label: "Several groups below", color: C.mixed, kind: "plain" });
    } else {
      legendTitle.textContent = `Tips coloured by ${META[f]} (n sequences)`;
      items = schemes[f].legend.concat([{ label: "No value", color: C.background, kind: "dot" }]);
    }
    items.push({ label: "No highlighted sequences", color: C.background, kind: "plain" });
    items.push({ label: "Collapsed group", color: C.group, kind: "patch" });
    legendEl.innerHTML = items.map(item =>
      `<span class="legend-item${item.g !== undefined ? " jump" : ""}"${item.g !== undefined ? ` data-g="${item.g}"` : ""}>${mark(item)}${esc(item.label)}</span>`).join("");
  }
  updateLegend();

  // ---------------------------------------------------------------
  // Layout and drawing
  // ---------------------------------------------------------------
  let matches = [], matchSet = new Set(), current = -1;
  const Y = new Float64Array(N.length);
  let geom = null;

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
        const h = c.length ? 1.6 : 1.0;
        Y[i] = row + h / 2; row += h; leaves.push(i);
        continue;
      }
      stack.push([i, true]);
      for (let k = c.length - 1; k >= 0; k--) stack.push([c[k], false]);
    }
    return { visible, leaves, rows: row };
  }

  function render() {
    if (!N.length) return;
    const L = layout();
    const W = wrap.clientWidth, rowH = +ui.row.value, zoom = +ui.zoom.value;
    const groupW = G.length ? (W >= 560 ? 130 : 100) : 0;
    const shown = [...ui.metaCols.querySelectorAll("input:checked")].map(el => +el.dataset.f);
    const metaW = W >= 560 ? 120 : 95;
    const colW = groupW + shown.length * metaW;
    const colorBy = +ui.colorBy.value;
    const labelRoom = Math.min(260, Math.max(140, W * 0.25));
    const fsHi = Math.max(9, Math.min(12, rowH * 0.8)), fsRef = fsHi - 1;

    // Horizontal window: from just before the highlighted MRCA to the furthest visible tip.
    const start = ui.truncate.checked ? X[D.mrca] : 0;
    let end = start;
    for (const i of L.leaves) end = Math.max(end, X[i]);
    const span = Math.max(end - start, 1e-9);
    const xStart = start - span * 0.04, triangleMax = span * 0.12;
    const treePx = Math.max(160, W - colW - labelRoom - 16) * zoom;
    const k = treePx / (end + triangleMax - xStart);
    const left = 8, top = 30;
    const svgW = Math.max(W, left + treePx + labelRoom + colW);
    const colX = svgW - colW + 6;
    const svgH = top + L.rows * rowH + 56;
    const px = x => left + (x - xStart) * k, py = y => top + y * rowH;
    const x0 = px(xStart);
    const bands = [], lines = [], shapes = [], texts = [], hits = [];

    for (const i of L.visible) {
      const n = N[i];
      const isLeaf = !n.c.length || collapsed.has(i);
      const x = isLeaf ? px(Math.max(X[i], xStart)) : px(X[i]);
      const y = py(Y[i]);

      if (!n.c.length) {
        if (matchSet.has(i)) {
          bands.push(`<rect x="0" y="${f1(y - rowH / 2)}" width="${svgW}" height="${rowH}" fill="${matches[current] === i ? "#ffd666" : "var(--match)"}"/>`);
        }
        const hi = n.g !== undefined;
        const mv = colorBy >= 0 ? metaValue(i, colorBy) : "";
        let fill = null;
        if (colorBy < 0) fill = hi ? GC[n.g] : null;
        else fill = mv ? schemes[colorBy].color(mv) : null;
        if (hi) {
          shapes.push(`<circle data-i="${i}" cx="${f1(x)}" cy="${f1(y)}" r="4" fill="${fill || C.background}" stroke="${colorBy < 0 ? C.surface : GC[n.g]}" stroke-width="${colorBy < 0 ? 1 : 2}"/>`);
          texts.push(`<text data-i="${i}" x="${f1(x + 7)}" y="${f1(y)}" font-size="${fsHi}" font-weight="700" fill="${C.text}" dominant-baseline="central">${esc(n.n)}</text>`);
          texts.push(`<text data-i="${i}" x="${f1(colX)}" y="${f1(y)}" font-size="${fsHi}" font-weight="700" fill="${GC[n.g]}" dominant-baseline="central">${esc(G[n.g])}</text>`);
        } else {
          shapes.push(`<circle data-i="${i}" cx="${f1(x)}" cy="${f1(y)}" ${fill ? `r="3" fill="${fill}" stroke="${C.surface}" stroke-width="0.8"` : `r="2" fill="${C.background}"`}/>`);
          texts.push(`<text data-i="${i}" x="${f1(x + 5)}" y="${f1(y)}" font-size="${fsRef}" fill="${C.muted}" dominant-baseline="central">${esc(n.n)}</text>`);
        }
        shown.forEach((f, slot) => {
          const v = metaValue(i, f);
          if (v) texts.push(`<text data-i="${i}" x="${f1(colX + groupW + slot * metaW)}" y="${f1(y)}" font-size="${hi ? fsHi : fsRef}" fill="${hi ? C.text2 : C.muted}" dominant-baseline="central">${esc(fit(v, metaW))}</text>`);
        });
        continue;
      }

      if (collapsed.has(i)) {
        const d = Math.min(Math.max(depth[i], span * 0.02), triangleMax) * k, h = 0.55 * rowH;
        shapes.push(`<g class="grp" data-i="${i}" data-act="expand"><polygon points="${f1(x)},${f1(y)} ${f1(x + d)},${f1(y - h)} ${f1(x + d)},${f1(y + h)}" fill="${C.group}" stroke="${C.background}" stroke-width="0.6"/>` +
          `<text x="${f1(x + d + 5)}" y="${f1(y)}" font-size="${fsRef}" fill="${C.muted}" dominant-baseline="central">${esc(collapsedLabel(i))}</text></g>`);
        continue;
      }

      const ys = n.c.map(c => py(Y[c]));
      lines.push(`<line x1="${f1(x)}" y1="${f1(Math.min(...ys))}" x2="${f1(x)}" y2="${f1(Math.max(...ys))}" stroke="${branchColor(n)}" stroke-width="${n.nh ? 1.4 : 1}"/>`);
      for (const c of n.c) {
        const yc = py(Y[c]);
        lines.push(`<line x1="${f1(x)}" y1="${f1(yc)}" x2="${f1(px(X[c]))}" y2="${f1(yc)}" stroke="${branchColor(N[c])}" stroke-width="${N[c].nh ? 1.8 : 1}"/>`);
        // Group name at the group's MRCA, when the group is a clade.
        if (N[c].l !== undefined) {
          const xc = collapsed.has(c) ? px(Math.max(X[c], xStart)) : px(X[c]);
          shapes.push(`<circle data-i="${c}" cx="${f1(xc)}" cy="${f1(yc)}" r="4.5" fill="${GC[N[c].l]}" stroke="${C.surface}" stroke-width="1"/>`);
          texts.push(`<text data-i="${c}" x="${f1(xc - 3)}" y="${f1(yc + rowH * 0.45)}" font-size="${fsHi}" font-weight="700" fill="${C.text}" text-anchor="end" dominant-baseline="hanging">${esc(G[N[c].l])}</text>`);
        }
      }
      if (i !== 0) hits.push(`<circle class="hit" data-i="${i}" data-act="collapse" cx="${f1(x)}" cy="${f1(y)}" r="6"/>`);
    }

    const length = Math.pow(10, Math.floor(Math.log10(span * 0.2)));
    const sbX = px(xStart + span * 0.02), sbY = svgH - 22;
    const scale = `<line x1="${f1(sbX)}" y1="${sbY}" x2="${f1(sbX + length * k)}" y2="${sbY}" stroke="${C.text2}" stroke-width="1.2"/>` +
      `<text x="${f1(sbX)}" y="${sbY - 5}" font-size="11" fill="${C.text2}">${+length.toPrecision(3)} subs/site</text>`;
    const head = (x, label) => `<text x="${f1(x)}" y="14" font-size="12" font-weight="700" fill="${C.text}" dominant-baseline="central">${esc(label)}</text>`;
    const headerTexts =
      (ui.truncate.checked ? `<text x="${left}" y="14" font-size="11" fill="${C.muted}" dominant-baseline="central">Root branch truncated at the MRCA of highlighted sequences</text>` : "") +
      (G.length ? head(colX, "Group") : "") +
      shown.map((f, slot) => head(colX + groupW + slot * metaW, cap(META[f]))).join("");

    wrap.innerHTML =
      `<svg id="tree" xmlns="http://www.w3.org/2000/svg" width="${f1(svgW)}" height="${f1(svgH)}" font-family='${FONT}'>` +
      `<rect width="100%" height="100%" fill="${C.surface}"/>` +
      `<defs><clipPath id="win"><rect x="${f1(x0)}" y="0" width="${f1(svgW - x0)}" height="${f1(svgH)}"/></clipPath></defs>` +
      bands.join("") +
      `<g clip-path="url(#win)" stroke-linecap="round">${lines.join("")}${shapes.join("")}${texts.join("")}${hits.join("")}</g>` +
      headerTexts + scale + `</svg>`;
    geom = { py };
  }

  // ---------------------------------------------------------------
  // Interaction
  // ---------------------------------------------------------------
  function scrollToNode(i) {
    if (geom) wrap.scrollTo({ top: Math.max(0, geom.py(Y[i]) - wrap.clientHeight / 2), behavior: "smooth" });
  }
  function reveal(i) { for (let j = parent[i]; j >= 0; j = parent[j]) collapsed.delete(j); }

  wrap.addEventListener("click", ev => {
    const el = ev.target.closest("[data-act]");
    if (!el) return;
    const i = +el.dataset.i;
    if (el.dataset.act === "expand") collapsed.delete(i); else collapsed.add(i);
    render();
  });

  legendEl.addEventListener("click", ev => {
    const el = ev.target.closest("[data-g]");
    if (!el) return;
    const g = +el.dataset.g;
    const first = N.findIndex(n => !n.c.length && n.g === g);
    if (first < 0) return;
    reveal(first); render(); scrollToNode(first);
  });

  const tip = document.getElementById("tip");
  wrap.addEventListener("mousemove", ev => {
    const el = ev.target.closest("[data-i]");
    if (!el) { tip.style.display = "none"; return; }
    const i = +el.dataset.i, n = N[i];
    let body;
    if (!n.c.length) {
      body = `<b>${esc(n.n)}</b>` + (n.g !== undefined ? `<br><span class="m">group</span> <b style="color:${GC[n.g]}">${esc(G[n.g])}</b>` : "") +
        META.map((name, f) => `<br><span class="m">${esc(name)}</span> ${esc(metaValue(i, f) || "—")}`).join("");
    } else {
      const groups = new Map();
      for (const t of tipsUnder(i)) if (N[t].g !== undefined) groups.set(N[t].g, (groups.get(N[t].g) || 0) + 1);
      body = (n.l !== undefined ? `<b>${esc(G[n.l])}</b><br>` : "") + `${n.nt} seqs · ${n.nh} highlighted` +
        (groups.size ? `<br><span class="m">groups</span> ${[...groups].map(([g, c]) => `${esc(G[g])} (${c})`).join(", ")}` : "") +
        (n.s ? `<br><span class="m">support</span> ${esc(n.s)}` : "") +
        META.map((label, f) => { const t = topValues(i, f, 3); return t ? `<br><span class="m">${esc(label)}</span> ${t}` : ""; }).join("") +
        (el.dataset.act || collapsed.has(i) ? `<br><span class="m">${collapsed.has(i) ? "Click to expand" : "Click to collapse"}</span>` : "");
    }
    tip.innerHTML = body;
    tip.style.display = "block";
    tip.style.left = Math.min(ev.clientX + 14, window.innerWidth - tip.offsetWidth - 8) + "px";
    tip.style.top = Math.min(ev.clientY + 14, window.innerHeight - tip.offsetHeight - 8) + "px";
  });
  wrap.addEventListener("mouseleave", () => { tip.style.display = "none"; });

  function runSearch() {
    const term = ui.search.value.trim().toLowerCase();
    matches = [];
    if (term.length >= 2) {
      for (let i = 0; i < N.length; i++) {
        const n = N[i];
        if (!n.c.length && (n.n.toLowerCase().includes(term) ||
            (n.g !== undefined && G[n.g].toLowerCase().includes(term)) ||
            (n.m || []).some(v => v.toLowerCase().includes(term)))) matches.push(i);
      }
    }
    matchSet = new Set(matches);
    matches.forEach(reveal);
    current = matches.length ? 0 : -1;
    ui.count.textContent = term.length < 2 ? "" : (matches.length ? `1 / ${matches.length}` : "0 found");
    render();
    if (current >= 0) scrollToNode(matches[current]);
  }
  let searchTimer = null;
  ui.search.addEventListener("input", () => { clearTimeout(searchTimer); searchTimer = setTimeout(runSearch, 200); });
  ui.search.addEventListener("keydown", ev => {
    if (ev.key !== "Enter" || !matches.length) return;
    ev.preventDefault();
    current = (current + (ev.shiftKey ? matches.length - 1 : 1)) % matches.length;
    ui.count.textContent = `${current + 1} / ${matches.length}`;
    render(); scrollToNode(matches[current]);
  });

  document.getElementById("expand-all").addEventListener("click", () => { collapsed = new Set(); render(); });
  document.getElementById("reset").addEventListener("click", () => { collapsed = defaultCollapsed(); matches.forEach(reveal); render(); });
  for (const el of [ui.truncate, ui.zoom, ui.row, ui.metaCols]) el.addEventListener("input", render);
  ui.colorBy.addEventListener("change", () => { updateLegend(); render(); });

  document.getElementById("download").addEventListener("click", () => {
    const svg = document.getElementById("tree");
    if (!svg) return;
    const blob = new Blob([new XMLSerializer().serializeToString(svg)], { type: "image/svg+xml" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = "highlight_tree.svg"; a.click();
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
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--tree", required=True, help="Newick tree file.")
    parser.add_argument("--metadata", required=True,
                        help="TSV/CSV with an accession column and the --columns below.")
    parser.add_argument("--highlight", required=True,
                        help="Accessions to highlight: ACCESSION<TAB>GROUP per line (GROUP optional).")
    parser.add_argument("--output", required=True, help="Output .html path.")
    parser.add_argument("--columns", nargs="+", default=["clade_name", "year", "host"],
                        help=("Metadata columns to show and colour by. A name also matches "
                              "a prefixed column (year -> collection_year). "
                              "Default: clade_name year host"))
    parser.add_argument("--id_column", default=None,
                        help=f"Accession column in --metadata. Default: first of {', '.join(ID_COLUMNS)}")
    parser.add_argument("--title", default=None, help="Page title. Default: tree file name.")
    return parser.parse_args()

TREE_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/rabv_tree.nwk"
METADATA_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/rabv_meta_data.csv"
OUTPUT_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/treeData.json"
def main():
    print("starting")
    # args = parse_args()
    # tree = Phylo.read(args.tree, "newick")
    tree = Phylo.read(TREE_FILE, "newick")
    tips = TipIndex([t.name for t in tree.get_terminals() if t.name])

    tip_group, groups = read_highlights(["clade_name", "year", "host"], tips)

    meta_columns, meta_values = read_metadata(METADATA_FILE, ["primary_accession", "EPA_major_clade"], "primary_accession", tips)
    nodes, mrca = tree_nodes(tree, tip_group, groups, meta_columns, meta_values)

    data = {
            "title": "test",
            "mrca": mrca,
            "groups": groups,
            "groupColors": group_colors(groups),
            "meta": meta_columns,
            "palette": PALETTE,
            "other": OTHER,
            "nodes": nodes,
        }

    with open(OUTPUT_FILE, "w", encoding="utf-8") as file:

      json.dump(
          data,
          file,
          indent=2
      )
    
    # title = args.title or Path(args.tree).name.split(".")[0]
    # path = write_html(args.output, title, nodes, mrca, groups, meta_columns)
    # print(f"Wrote {path}  ({len(tips.exact)} tips)")
    return 0


if __name__ == "__main__":
  main()
#     # Large ladder-like trees need deep recursion in Bio.Phylo (common_ancestor, get_terminals).
#     sys.setrecursionlimit(1_000_000)
#     threading.stack_size(512 * 1024 * 1024)
#     result = []
#     worker = threading.Thread(target=lambda: result.append(main()))
#     worker.start()
#     worker.join()
#     raise SystemExit(result[0] if result else 1)
