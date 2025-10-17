# import sqlite3
# import csv

# # 1. Connect to database
# conn = sqlite3.connect("./db/rabv-gdb-jul092025.db")
# cur = conn.cursor()

# # 2. Add new columns (if not already there)
# try:
#     cur.execute("ALTER TABLE meta_data ADD COLUMN major_clade TEXT;")
# except sqlite3.OperationalError:
#     pass  # column already exists

# try:
#     cur.execute("ALTER TABLE meta_data ADD COLUMN minor_clade TEXT;")
# except sqlite3.OperationalError:
#     pass  # column already exists

# # 3. Read CSV and update rows
# with open("/Users/danaallen/Downloads/for_dana/major_clades.tsv", newline="") as csvfile:
#     reader = csv.DictReader(csvfile, delimiter="\t")
#     for row in reader:
#         cur.execute("""
#             UPDATE meta_data
#             SET major_clade = ?
#             WHERE primary_accession = ?;
#         """, (row["clade"], row["acc"]))

#         print(f"Updating acc={row["acc"]}, clade={row["clade"]}, rows affected={cur.rowcount}")

# # 3. Read CSV and update rows
# with open("/Users/danaallen/Downloads/for_dana/minor_clades.tsv", newline="") as csvfile:
#     reader = csv.DictReader(csvfile, delimiter="\t")
#     for row in reader:
#         cur.execute("""
#             UPDATE meta_data
#             SET minor_clade = ?
#             WHERE primary_accession = ?;
#         """, (row["subclade"], row["acc"]))

# # 4. Save changes and close
# conn.commit()
# conn.close()


import sqlite3

conn = sqlite3.connect("./db/rabv-gdb-jul092025.db")
cur = conn.cursor()

# Drop old table if it exists
cur.execute("DROP TABLE IF EXISTS genotypes;")

# Create new genotypes table
cur.execute("""
    CREATE TABLE genotypes AS
    -- case 1: majors with minors
    SELECT DISTINCT major_clade,
                    minor_clade
    FROM meta_data
    WHERE major_clade IS NOT NULL
      AND minor_clade IS NOT NULL

    UNION

    -- case 2: majors with no minors
    SELECT DISTINCT major_clade,
                    NULL AS minor_clade
    FROM meta_data
    WHERE major_clade IS NOT NULL
      AND minor_clade IS NULL;
""")

conn.commit()
conn.close()