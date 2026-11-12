
# import json
# import csv
# from Bio import Phylo

# # --------------------------------------------------
# # INPUT FILES
# # --------------------------------------------------

TREE_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/rabv_tree.nwk"
METADATA_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/rabv_meta_data.csv"
OUTPUT_FILE = "/Users/dana/CVR/gDB/backend/V-gDB/data_import/hcv/treeData.json"

# ACCESSION_COLUMN = "primary_accession"

import csv
import json
from Bio import Phylo

# --------------------------------------------------
# CONFIGURATION
# --------------------------------------------------

ACCESSION_COLUMN = "primary_accession"
MAJOR_COLUMN = "EPA_major_clade"
MINOR_COLUMN = "EPA_minor_clade"

META_COLUMNS = [
    "country",
    "host",
    "collection_year",
]

# --------------------------------------------------
# MAJOR AND MINOR CLADE DEFINITIONS
# --------------------------------------------------

# Use the same order as your existing JSON arrays.
# Paste your complete existing arrays here.
majors = [
    {"name": "Africa-2", "color": "#008300", "n": 892},
    {"name": "Africa-3", "color": "#4a3aa7", "n": 368},
    {"name": "Arctic", "color": "#eda100", "n": 3225},
    {"name": "Asian", "color": "#1baf7a", "n": 5352},
    {"name": "Bats", "color": "#eb6834", "n": 7515},
    {"name": "Cosmopolitan", "color": "#2a78d6", "n": 12321},
    {"name": "Indian-Sub", "color": "#e34948", "n": 314},
    {"name": "RAC-SK", "color": "#e87ba4", "n": 2315},
    {"name": "Unclassified", "color": "#8a8984", "n": 82},
]

# IMPORTANT:
# Insert your complete subclades array here, in its original order.
subclades = [
    {
      "name": "A",
      "major": 2,
      "n": 921
    },
    {
      "name": "AF1a",
      "major": 5,
      "n": 1874
    },
    {
      "name": "AF1b",
      "major": 5,
      "n": 3593
    },
    {
      "name": "AF1c",
      "major": 5,
      "n": 40
    },
    {
      "name": "AF4",
      "major": 5,
      "n": 53
    },
    {
      "name": "AL1a",
      "major": 2,
      "n": 1823
    },
    {
      "name": "AL1b",
      "major": 2,
      "n": 132
    },
    {
      "name": "AL2",
      "major": 2,
      "n": 257
    },
    {
      "name": "AL3",
      "major": 2,
      "n": 91
    },
    {
      "name": "AM1",
      "major": 5,
      "n": 354
    },
    {
      "name": "AM2a",
      "major": 5,
      "n": 773
    },
    {
      "name": "AM2b",
      "major": 5,
      "n": 55
    },
    {
      "name": "AM3a",
      "major": 5,
      "n": 354
    },
    {
      "name": "AM3b",
      "major": 5,
      "n": 218
    },
    {
      "name": "AM3c",
      "major": 5,
      "n": 61
    },
    {
      "name": "AM4",
      "major": 5,
      "n": 38
    },
    {
      "name": "AM5",
      "major": 5,
      "n": 193
    },
    {
      "name": "AP",
      "major": 4,
      "n": 12
    },
    {
      "name": "Africa-2",
      "major": 0,
      "n": 892
    },
    {
      "name": "Africa-3",
      "major": 1,
      "n": 368
    },
    {
      "name": "Asian",
      "major": 3,
      "n": 10
    },
    {
      "name": "Bats",
      "major": 4,
      "n": 19
    },
    {
      "name": "CA1",
      "major": 5,
      "n": 935
    },
    {
      "name": "CA2",
      "major": 5,
      "n": 361
    },
    {
      "name": "CA3",
      "major": 5,
      "n": 207
    },
    {
      "name": "CE",
      "major": 5,
      "n": 178
    },
    {
      "name": "Cosmopolitan",
      "major": 5,
      "n": 4
    },
    {
      "name": "DR",
      "major": 4,
      "n": 4620
    },
    {
      "name": "EE",
      "major": 5,
      "n": 234
    },
    {
      "name": "EF-E1",
      "major": 4,
      "n": 119
    },
    {
      "name": "EF-E2",
      "major": 4,
      "n": 449
    },
    {
      "name": "EF-W1",
      "major": 4,
      "n": 227
    },
    {
      "name": "EF-W2",
      "major": 4,
      "n": 80
    },
    {
      "name": "EFu",
      "major": 4,
      "n": 99
    },
    {
      "name": "Indian-Sub",
      "major": 6,
      "n": 314
    },
    {
      "name": "LB1",
      "major": 4,
      "n": 30
    },
    {
      "name": "LB2",
      "major": 4,
      "n": 106
    },
    {
      "name": "LC",
      "major": 4,
      "n": 237
    },
    {
      "name": "LI",
      "major": 4,
      "n": 102
    },
    {
      "name": "LN",
      "major": 4,
      "n": 128
    },
    {
      "name": "LS",
      "major": 4,
      "n": 22
    },
    {
      "name": "LX",
      "major": 4,
      "n": 19
    },
    {
      "name": "ME1a",
      "major": 5,
      "n": 565
    },
    {
      "name": "ME1b",
      "major": 5,
      "n": 17
    },
    {
      "name": "ME2",
      "major": 5,
      "n": 275
    },
    {
      "name": "MYsp",
      "major": 4,
      "n": 203
    },
    {
      "name": "MYu",
      "major": 4,
      "n": 355
    },
    {
      "name": "NEE",
      "major": 5,
      "n": 576
    },
    {
      "name": "PH",
      "major": 4,
      "n": 99
    },
    {
      "name": "PS",
      "major": 4,
      "n": 51
    },
    {
      "name": "RAC-SK",
      "major": 7,
      "n": 2315
    },
    {
      "name": "SEA1a",
      "major": 3,
      "n": 37
    },
    {
      "name": "SEA1b",
      "major": 3,
      "n": 2144
    },
    {
      "name": "SEA2a",
      "major": 3,
      "n": 337
    },
    {
      "name": "SEA2b",
      "major": 3,
      "n": 74
    },
    {
      "name": "SEA3",
      "major": 3,
      "n": 1275
    },
    {
      "name": "SEA4",
      "major": 3,
      "n": 1323
    },
    {
      "name": "SEA5",
      "major": 3,
      "n": 152
    },
    {
      "name": "TB1",
      "major": 4,
      "n": 157
    },
    {
      "name": "TB2",
      "major": 4,
      "n": 381
    },
    {
      "name": "Unclassified",
      "major": 8,
      "n": 82
    },
    {
      "name": "Vac",
      "major": 5,
      "n": 543
    },
    {
      "name": "Vac2",
      "major": 5,
      "n": 313
    },
    {
      "name": "WE",
      "major": 5,
      "n": 448
    },
    {
      "name": "YUGCOW",
      "major": 5,
      "n": 8
    },
    {
      "name": "YUGFOX",
      "major": 5,
      "n": 52
    }
  ]

# ==================================================
# BUILD LOOKUP DICTIONARIES
# ==================================================

major_lookup = {
    item["name"]: index
    for index, item in enumerate(majors)
}

subclade_lookup = {
    item["name"]: index
    for index, item in enumerate(subclades)
}

UNCLASSIFIED_MAJOR = major_lookup["Unclassified"]


# ==================================================
# READ TREE
# ==================================================

tree = Phylo.read(TREE_FILE, "newick")


# ==================================================
# READ METADATA
# ==================================================

with open(METADATA_FILE, newline="", encoding="utf-8-sig") as file:

    reader = csv.DictReader(file)

    metadata = {
        row[ACCESSION_COLUMN].strip(): row
        for row in reader
        if row.get(ACCESSION_COLUMN)
        and row[ACCESSION_COLUMN].strip()
    }


# ==================================================
# BUILD METADATA VALUE LOOKUPS
# ==================================================

meta_values = []
meta_lookup = []

for column in META_COLUMNS:

    values = sorted({
        row[column].strip()
        for row in metadata.values()
        if row.get(column) and row[column].strip()
    })

    values = [""] + values

    meta_values.append(values)

    meta_lookup.append({
        value: index
        for index, value in enumerate(values)
    })


# ==================================================
# CONVERT TREE TO COMPACT NODE FORMAT
# ==================================================

nodes = []


def add_node(clade):

    node_index = len(nodes)
    nodes.append(None)

    # --------------------------------------------------
    # RECURSIVELY BUILD CHILDREN
    # --------------------------------------------------

    children = []

    for child_clade in clade.clades:

        child_index = add_node(child_clade)
        children.append(child_index)

    # --------------------------------------------------
    # SEQUENCE COUNT
    # --------------------------------------------------

    if clade.is_terminal():

        sequence_count = 1

    else:

        sequence_count = sum(
            nodes[index]["nt"]
            for index in children
        )

    # --------------------------------------------------
    # MAJOR AND MINOR CLADE ASSIGNMENTS
    # --------------------------------------------------

    major_index = -1
    subclade_index = -1

    if clade.is_terminal():

        accession = clade.name
        row = metadata.get(accession, {})

        major_name = (
            row.get(MAJOR_COLUMN) or ""
        ).strip()

        minor_name = (
            row.get(MINOR_COLUMN) or ""
        ).strip()

        # Major clade
        major_index = major_lookup.get(
            major_name,
            UNCLASSIFIED_MAJOR
        )

        # Minor clade
        candidate_subclade = subclade_lookup.get(
            minor_name,
            -1
        )

        if candidate_subclade != -1:

            expected_major = subclades[candidate_subclade]["major"]

            if expected_major == major_index:
                subclade_index = candidate_subclade

    else:

        # Determine whether all children belong
        # to the same major/minor clade.

        child_majors = [
            nodes[index]["mj"]
            for index in children
        ]

        child_subclades = [
            nodes[index]["sc"]
            for index in children
        ]

        # Major clade inheritance
        if (
            child_majors
            and child_majors[0] != -1
            and all(
                value == child_majors[0]
                for value in child_majors
            )
        ):

            major_index = child_majors[0]

        # Minor clade inheritance
        if (
            child_subclades
            and child_subclades[0] != -1
            and all(
                value == child_subclades[0]
                for value in child_subclades
            )
        ):

            subclade_index = child_subclades[0]

    # --------------------------------------------------
    # BUILD NODE
    # --------------------------------------------------

    node = {
        "c": children,
        "b": clade.branch_length or 0.0,
        "nt": sequence_count,
        "mj": major_index,
        "sc": subclade_index,
    }

    # Branch support
    if clade.confidence is not None:
        node["s"] = str(clade.confidence)

    # --------------------------------------------------
    # TERMINAL METADATA
    # --------------------------------------------------

    if clade.is_terminal():

        accession = clade.name
        node["n"] = accession

        row = metadata.get(accession, {})

        node["m"] = [
            meta_lookup[i].get(
                row.get(column, "").strip(),
                0
            )
            for i, column in enumerate(META_COLUMNS)
        ]

    nodes[node_index] = node

    return node_index


# ==================================================
# RUN CONVERSION
# ==================================================

add_node(tree.root)


# ==================================================
# BUILD PARENT LOOKUP
# ==================================================

parent = [-1] * len(nodes)

for i, node in enumerate(nodes):

    for child in node["c"]:
        parent[child] = i


# ==================================================
# IDENTIFY MAJOR CLADE ROOTS (pm)
# ==================================================

major_clade_roots = []

for i, node in enumerate(nodes):

    parent_index = parent[i]

    # Skip the root
    if parent_index == -1:
        continue

    # Only internal nodes can represent a major-clade root
    if not node["c"]:
        continue

    parent_major = nodes[parent_index]["mj"]
    current_major = node["mj"]

    # Mark the beginning of a new classified major clade
    if (
        current_major >= 0
        and current_major != parent_major
    ):

        node["pm"] = 1
        major_clade_roots.append(i)

# ==================================================
# IDENTIFY SUBCLADE ROOTS (pt)
# ==================================================

subclade_roots = []

for i, node in enumerate(nodes):

    parent_index = parent[i]

    # Skip the root
    if parent_index == -1:
        continue

    # Only internal nodes can represent a subclade root
    if not node["c"]:
        continue

    parent_subclade = nodes[parent_index]["sc"]
    current_subclade = node["sc"]

    # Mark the beginning of a new classified subclade
    if (
        current_subclade >= 0
        and current_subclade != parent_subclade
    ):

        node["pt"] = 1
        subclade_roots.append(i)

# ==================================================
# BUILD OUTPUT
# ==================================================

output = {
    "meta": META_COLUMNS,
    "metaValues": meta_values,
    "majors": majors,
    "subclades": subclades,
    "nodes": nodes,
}


# ==================================================
# EXPORT JSON
# ==================================================

with open(OUTPUT_FILE, "w", encoding="utf-8") as file:

    json.dump(
        output,
        file,
        indent=2
    )


# ==================================================
# SUMMARY
# ==================================================

print("Tree exported successfully")

print(f"Nodes: {len(nodes):,}")
print(f"Sequences: {nodes[0]['nt']:,}")
print(f"Major clades: {len(majors)}")
print(f"Subclades: {len(subclades)}")
print(f"Major clade roots (pm): {len(major_clade_roots):,}")


# ==================================================
# CHECK TREE CONNECTIVITY
# ==================================================

def check_tree(nodes):

    visited = set()
    stack = [0]

    parent_count = [0] * len(nodes)

    for node in nodes:

        for child in node.get("c", []):

            if child < 0 or child >= len(nodes):
                print("Invalid child index:", child)
                continue

            parent_count[child] += 1

    while stack:

        index = stack.pop()

        if index in visited:
            continue

        visited.add(index)

        stack.extend(nodes[index].get("c", []))

    print("\nTree connectivity")

    print("Total nodes:", len(nodes))
    print("Reachable from root:", len(visited))
    print("Disconnected nodes:", len(nodes) - len(visited))

    print(
        "Nodes with multiple parents:",
        [
            i for i, count in enumerate(parent_count)
            if count > 1
        ]
    )

    print(
        "Nodes without parents (excluding root):",
        [
            i for i, count in enumerate(parent_count)
            if count == 0 and i != 0
        ][:20]
    )

    print(
        "Nodes marked pm:",
        sum(1 for node in nodes if node.get("pm") == 1)
    )


check_tree(nodes)

# # Metadata fields to include in the tree
# META_COLUMNS = [
#     "country",
#     "host",
#     "collection_year",
# ]
# import csv
# from Bio import Phylo

# # --------------------------------------------------
# # CONFIGURATION
# # --------------------------------------------------
# ACCESSION_COLUMN = "primary_accession"

# META_COLUMNS = [
#     "country",
#     "host",
#     "collection_year",
# ]

# # --------------------------------------------------
# # READ TREE
# # --------------------------------------------------

# tree = Phylo.read(TREE_FILE, "newick")

# # --------------------------------------------------
# # READ METADATA
# # --------------------------------------------------

# with open(METADATA_FILE, newline="", encoding="utf-8-sig") as file:
#     reader = csv.DictReader(file)

#     metadata = {
#         row[ACCESSION_COLUMN].strip(): row
#         for row in reader
#         if row.get(ACCESSION_COLUMN)
#     }

# # --------------------------------------------------
# # BUILD METADATA VALUE LOOKUPS
# # --------------------------------------------------

# # Use 1-based indexes, reserving index 0.
# meta_values = []
# meta_lookup = []

# for column in META_COLUMNS:

#     values = sorted({
#         row[column].strip()
#         for row in metadata.values()
#         if row.get(column) and row[column].strip()
#     })

#     values = [""] + values

#     meta_values.append(values)

#     meta_lookup.append({
#         value: index
#         for index, value in enumerate(values)
#     })

# # --------------------------------------------------
# # BUILD TREE ANNOTATIONS
# # --------------------------------------------------

# def annotate_tree(clade):

#     # Number of terminal descendants
#     if clade.is_terminal():
#         clade.nt = 1
#     else:
#         clade.nt = sum(
#             annotate_tree(child)
#             for child in clade.clades
#         )

#     # Preserve branch support when available
#     clade.s = (
#         str(clade.confidence)
#         if clade.confidence is not None
#         else None
#     )

#     # Default classification values
#     clade.mj = -1
#     clade.sc = -1

#     # Additional annotation placeholders
#     clade.pm = None
#     clade.pt = None
#     clade.dp = None
#     clade.ps = None

#     return clade.nt


# annotate_tree(tree.root)

# # --------------------------------------------------
# # CONVERT TREE TO COMPACT NODE FORMAT
# # --------------------------------------------------

# nodes = []


# def add_node(clade):

#     node_index = len(nodes)

#     # Reserve position before processing children.
#     # This preserves preorder indexing.
#     nodes.append(None)

#     children = []

#     for child_clade in clade.clades:
#         child_index = add_node(child_clade)
#         children.append(child_index)

#     node = {
#         "c": children,
#         "b": clade.branch_length or 0.0,
#         "nt": clade.nt,
#         "mj": clade.mj,
#         "sc": clade.sc,
#     }

#     # Branch support
#     if clade.s is not None:
#         node["s"] = clade.s

#     # Additional annotations, when available
#     for field in ("pm", "pt", "dp", "ps"):
#         value = getattr(clade, field, None)

#         if value is not None:
#             node[field] = value

#     # Terminal sequence metadata
#     if clade.is_terminal():

#         accession = clade.name
#         node["n"] = accession

#         row = metadata.get(accession, {})

#         node["m"] = [
#             meta_lookup[i].get(
#                 row.get(column, "").strip(),
#                 0
#             )
#             for i, column in enumerate(META_COLUMNS)
#         ]

#     nodes[node_index] = node

#     return node_index


# add_node(tree.root)

# # --------------------------------------------------
# # BUILD FINAL JSON
# # --------------------------------------------------

# import json

# # output = {
# #     "meta": META_COLUMNS,
# #     "metaValues": meta_values,
# #     "nodes": nodes,
# # }

# # with open(OUTPUT_FILE, "w", encoding="utf-8") as file:
# #     json.dump(output, file, indent=2)

# print(f"Tree exported successfully: {OUTPUT_FILE}")
# print(f"Total nodes: {len(nodes):,}")
# print(f"Total sequences: {tree.root.nt:,}")

# tree_data = {
#     "title": "Rabies Virus Phylogenetic Tree",
#     "banner": "",
#     "purity": 0,
#     "unknown": 0,
#     "nLineages": 0,
#     "meta": META_COLUMNS,
#     "metaValues": meta_values,
#     "lin": [""],
#     "subclades": [
#     {
#       "name": "A",
#       "major": 2,
#       "n": 921
#     },
#     {
#       "name": "AF1a",
#       "major": 5,
#       "n": 1874
#     },
#     {
#       "name": "AF1b",
#       "major": 5,
#       "n": 3593
#     },
#     {
#       "name": "AF1c",
#       "major": 5,
#       "n": 40
#     },
#     {
#       "name": "AF4",
#       "major": 5,
#       "n": 53
#     },
#     {
#       "name": "AL1a",
#       "major": 2,
#       "n": 1823
#     },
#     {
#       "name": "AL1b",
#       "major": 2,
#       "n": 132
#     },
#     {
#       "name": "AL2",
#       "major": 2,
#       "n": 257
#     },
#     {
#       "name": "AL3",
#       "major": 2,
#       "n": 91
#     },
#     {
#       "name": "AM1",
#       "major": 5,
#       "n": 354
#     },
#     {
#       "name": "AM2a",
#       "major": 5,
#       "n": 773
#     },
#     {
#       "name": "AM2b",
#       "major": 5,
#       "n": 55
#     },
#     {
#       "name": "AM3a",
#       "major": 5,
#       "n": 354
#     },
#     {
#       "name": "AM3b",
#       "major": 5,
#       "n": 218
#     },
#     {
#       "name": "AM3c",
#       "major": 5,
#       "n": 61
#     },
#     {
#       "name": "AM4",
#       "major": 5,
#       "n": 38
#     },
#     {
#       "name": "AM5",
#       "major": 5,
#       "n": 193
#     },
#     {
#       "name": "AP",
#       "major": 4,
#       "n": 12
#     },
#     {
#       "name": "Africa-2",
#       "major": 0,
#       "n": 892
#     },
#     {
#       "name": "Africa-3",
#       "major": 1,
#       "n": 368
#     },
#     {
#       "name": "Asian",
#       "major": 3,
#       "n": 10
#     },
#     {
#       "name": "Bats",
#       "major": 4,
#       "n": 19
#     },
#     {
#       "name": "CA1",
#       "major": 5,
#       "n": 935
#     },
#     {
#       "name": "CA2",
#       "major": 5,
#       "n": 361
#     },
#     {
#       "name": "CA3",
#       "major": 5,
#       "n": 207
#     },
#     {
#       "name": "CE",
#       "major": 5,
#       "n": 178
#     },
#     {
#       "name": "Cosmopolitan",
#       "major": 5,
#       "n": 4
#     },
#     {
#       "name": "DR",
#       "major": 4,
#       "n": 4620
#     },
#     {
#       "name": "EE",
#       "major": 5,
#       "n": 234
#     },
#     {
#       "name": "EF-E1",
#       "major": 4,
#       "n": 119
#     },
#     {
#       "name": "EF-E2",
#       "major": 4,
#       "n": 449
#     },
#     {
#       "name": "EF-W1",
#       "major": 4,
#       "n": 227
#     },
#     {
#       "name": "EF-W2",
#       "major": 4,
#       "n": 80
#     },
#     {
#       "name": "EFu",
#       "major": 4,
#       "n": 99
#     },
#     {
#       "name": "Indian-Sub",
#       "major": 6,
#       "n": 314
#     },
#     {
#       "name": "LB1",
#       "major": 4,
#       "n": 30
#     },
#     {
#       "name": "LB2",
#       "major": 4,
#       "n": 106
#     },
#     {
#       "name": "LC",
#       "major": 4,
#       "n": 237
#     },
#     {
#       "name": "LI",
#       "major": 4,
#       "n": 102
#     },
#     {
#       "name": "LN",
#       "major": 4,
#       "n": 128
#     },
#     {
#       "name": "LS",
#       "major": 4,
#       "n": 22
#     },
#     {
#       "name": "LX",
#       "major": 4,
#       "n": 19
#     },
#     {
#       "name": "ME1a",
#       "major": 5,
#       "n": 565
#     },
#     {
#       "name": "ME1b",
#       "major": 5,
#       "n": 17
#     },
#     {
#       "name": "ME2",
#       "major": 5,
#       "n": 275
#     },
#     {
#       "name": "MYsp",
#       "major": 4,
#       "n": 203
#     },
#     {
#       "name": "MYu",
#       "major": 4,
#       "n": 355
#     },
#     {
#       "name": "NEE",
#       "major": 5,
#       "n": 576
#     },
#     {
#       "name": "PH",
#       "major": 4,
#       "n": 99
#     },
#     {
#       "name": "PS",
#       "major": 4,
#       "n": 51
#     },
#     {
#       "name": "RAC-SK",
#       "major": 7,
#       "n": 2315
#     },
#     {
#       "name": "SEA1a",
#       "major": 3,
#       "n": 37
#     },
#     {
#       "name": "SEA1b",
#       "major": 3,
#       "n": 2144
#     },
#     {
#       "name": "SEA2a",
#       "major": 3,
#       "n": 337
#     },
#     {
#       "name": "SEA2b",
#       "major": 3,
#       "n": 74
#     },
#     {
#       "name": "SEA3",
#       "major": 3,
#       "n": 1275
#     },
#     {
#       "name": "SEA4",
#       "major": 3,
#       "n": 1323
#     },
#     {
#       "name": "SEA5",
#       "major": 3,
#       "n": 152
#     },
#     {
#       "name": "TB1",
#       "major": 4,
#       "n": 157
#     },
#     {
#       "name": "TB2",
#       "major": 4,
#       "n": 381
#     },
#     {
#       "name": "Unclassified",
#       "major": 8,
#       "n": 82
#     },
#     {
#       "name": "Vac",
#       "major": 5,
#       "n": 543
#     },
#     {
#       "name": "Vac2",
#       "major": 5,
#       "n": 313
#     },
#     {
#       "name": "WE",
#       "major": 5,
#       "n": 448
#     },
#     {
#       "name": "YUGCOW",
#       "major": 5,
#       "n": 8
#     },
#     {
#       "name": "YUGFOX",
#       "major": 5,
#       "n": 52
#     }
#   ],
#   "majors": [
#     {
#       "name": "Africa-2",
#       "color": "#008300",
#       "n": 892
#     },
#     {
#       "name": "Africa-3",
#       "color": "#4a3aa7",
#       "n": 368
#     },
#     {
#       "name": "Arctic",
#       "color": "#eda100",
#       "n": 3225
#     },
#     {
#       "name": "Asian",
#       "color": "#1baf7a",
#       "n": 5352
#     },
#     {
#       "name": "Bats",
#       "color": "#eb6834",
#       "n": 7515
#     },
#     {
#       "name": "Cosmopolitan",
#       "color": "#2a78d6",
#       "n": 12321
#     },
#     {
#       "name": "Indian-Sub",
#       "color": "#e34948",
#       "n": 314
#     },
#     {
#       "name": "RAC-SK",
#       "color": "#e87ba4",
#       "n": 2315
#     },
#     {
#       "name": "Unclassified",
#       "color": "#8a8984",
#       "n": 82
#     }
#   ],
#     "nodes": nodes,
# }

# with open(OUTPUT_FILE, "w", encoding="utf-8") as file:
#     json.dump(tree_data, file, separators=(",", ":"))

# # print(f"Generated {OUTPUT_FILE}")
# # print(f"Total nodes: {len(nodes)}")
# # print(f"Sequences: {nodes[0]['nt']}")
# # print(f"Metadata matches: {sum(1 for node in nodes if node and node.get('n') in metadata)}")
