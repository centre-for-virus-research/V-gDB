# import sqlite3
import csv

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

conn = sqlite3.connect("/Volumes/My Passport/V-gDB-flu-Sep292025.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

# # Drop old table if it exists
cur.execute("DROP TABLE IF EXISTS adaptive_mutations;")

# Create new genotypes table
cur.execute("""CREATE TABLE adaptive_mutations 
            (
              segment TEXT,
              mutation TEXT,
              wt TEXT,
              position TEXT,
              mutant TEXT,
              authors TEXT,
              experimentally_verified TEXT,
              discovery TEXT,
              virus TEXT,
              notes TEXT,
              doi TEXT
            )""")

conn.commit()

with open("/Volumes/My Passport/CVR/gdb/Flu/adaptive_mutations/Clean M.csv", newline="") as csvfile:
    reader = csv.DictReader(csvfile, delimiter=",")
    for row in reader:
        cur.execute("""
            INSERT INTO adaptive_mutations (
              segment,
              mutation,
              wt,
              position,
              mutant,
              authors,
              experimentally_verified,
              discovery,
              virus,
              notes, 
              doi
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            row['Segment'],
            row['Mutation'],
            row['WT'],
            row['Position'],
            row['Mutant'],
            row['Paper..s.'],
            row['Experimentally.Verified'],
            row['How.was.the.mutation.discovered'],
            row['Which.virus.'],
            row['Notes'],
            row['DOI']
        ))

#         print(f"Updating acc={row["acc"]}, clade={row["clade"]}, rows affected={cur.rowcount}")
  
  
conn.commit()
conn.close()


