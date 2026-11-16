from __future__ import annotations

import argparse
import csv
import html
import json
import sqlite3
import sys
from pathlib import Path
# import phylotreelib as pt
from Bio import Phylo
from io import StringIO
import re
from django.db import connections
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
# VERSION = re.compile(r"\.\d+$")
# class TipIndex:
#     def __init__(self, names):
#         self.exact = set(names)
#         self.base = {}
#         for name in names:
#             self.base.setdefault(VERSION.sub("", name), name)

#     def find(self, accession):
#         accession = accession.strip()
#         if accession in self.exact:
#             return accession
#         return self.base.get(VERSION.sub("", accession))



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


def read_clades(tips, db):
    """{tip: (major, minor)}: EPA tables of the references first, then the database."""
    out = {}
    if db:
      with connections[db].cursor() as cursor:

          cursor.execute("""
              SELECT
                  primary_accession,
                  accession_version,
                  EPA_major_clade,
                  EPA_minor_clade
              FROM meta_data
          """)

          for acc, version, ma, mi in cursor:
              tip = tips.find(acc or "") or tips.find(version or "")

              if tip and tip not in out and clean(ma):
                  out[tip] = (clean(ma), clean(mi))
    return out


# def read_lineages(lineage_dir, tip_sub, tips):
#     """{tip: (lineage, parent)} from each sub-clade's lineage table; a sequence only takes
#     the lineage from the table of its own sub-clade."""
#     out, missing = {}, []
#     for sub in sorted(set(tip_sub.values())):
#         if sub == UNCLASSIFIED:
#             continue
#         base = Path(lineage_dir) / sub
#         path = next((p for p in (base / "meta_data.tsv", base / f"{sub}.tip_lineage_assignments.tsv")
#                      if p.exists()), None)
#         if path is None:
#             missing.append(sub)
#             continue
#         with path.open(newline="", encoding="utf-8") as handle:
#             for row in csv.DictReader(handle, delimiter="\t"):
#                 tip = tips.find(row.get("ID") or row.get("sequence_id") or "")
#                 if tip and tip_sub.get(tip) == sub:
#                     out[tip] = (clean(row.get("Lineage") or row.get("lineage")),
#                                 clean(row.get("Parent") or row.get("lineage_parent")))
#     if missing:
#         print(f"Warning: no lineage table in {lineage_dir} for: {', '.join(missing)}")
#     return out
def read_lineages(db, tip_sub, tips):
    """Return {tip: (lineage, (major_clade, minor_clade))}."""

    out = {}

    with connections[db].cursor() as cursor:
        cursor.execute("""
            SELECT
                primary_accession,
                EPA_major_clade,
                EPA_minor_clade,
                maddog_lineage
            FROM meta_data
            WHERE maddog_lineage IS NOT NULL
              AND maddog_lineage != ''
        """)

        for acc, major, minor, lineage in cursor:
            tip = tips.find(acc or "")

            if not tip:
                continue

            if tip_sub.get(tip) != clean(minor):
                continue

            out[tip] = (
                clean(lineage),
                (clean(major), clean(minor))
            )

    return out

# def read_db_metadata(db, tips):
#     con = sqlite3.connect(db)
#     values = {c: {} for c in META_COLUMNS}
#     query = "SELECT primary_accession, accession_version, {} FROM meta_data".format(", ".join(META_COLUMNS))
#     for acc, version, *row in con.execute(query):
#         tip = tips.find(acc or "") or tips.find(version or "")
#         if tip is None:
#             continue
#         for column, value in zip(META_COLUMNS, row):
#             value = clean(value)
#             if value:
#                 values[column].setdefault(tip, value)
#     con.close()
#     return values
def read_db_metadata(db, tips):
    values = {c: {} for c in META_COLUMNS}

    query = """
        SELECT
            primary_accession,
            accession_version,
            {}
        FROM meta_data
    """.format(", ".join(META_COLUMNS))

    with connections[db].cursor() as cursor:
        cursor.execute(query)

        for acc, version, *row in cursor:
            tip = tips.find(acc or "") or tips.find(version or "")

            if tip is None:
                continue

            for column, value in zip(META_COLUMNS, row):
                value = clean(value)

                if value:
                    values[column].setdefault(tip, value)

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

def build_maddog_tree(newick, db):
    
    tree = Phylo.read(StringIO(newick), "newick")
    tips = TipIndex([t.name for t in tree.get_terminals() if t.name])

    clades = read_clades(tips, db)
    tip_major, tip_sub = {}, {}
    for tip in tips.exact:
        major, minor = clades.get(tip, ("", ""))
        major = major or UNCLASSIFIED
        tip_major[tip] = major
        tip_sub[tip] = minor if minor and major != UNCLASSIFIED else major
    unclassified = sum(1 for v in tip_major.values() if v == UNCLASSIFIED)
    if unclassified:
        print(f"Note: {unclassified} tips have no EPA major clade; shown as '{UNCLASSIFIED}'")

    lineages = read_lineages(db, tip_sub, tips)
    meta_values = read_db_metadata(db, tips)

    nodes, end = flat_nodes(tree)
    subclades, majors = annotate(nodes, end, tip_major, tip_sub, lineages, 0.9)
    lin_list, meta_lists = encode_tips(nodes, lineages, meta_values)
    n_lineages = sum(1 for v in lineages.values() if v[0] and v[0] != "Unassigned")
    data = {
        "title": "test", "banner": BANNER, "purity": 0.9, "unknown": UNKNOWN,
        "nLineages": n_lineages, "meta": META_COLUMNS, "metaValues": meta_lists, "lin": lin_list,
        "subclades": subclades, "majors": majors, "nodes": nodes,
    }

    return data
