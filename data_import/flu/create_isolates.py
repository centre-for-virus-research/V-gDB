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

conn = sqlite3.connect("/Volumes/My Passport/V-gDB-flu-Sep292025.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

# # Drop old table if it exists
cur.execute("DROP TABLE IF EXISTS isolates;")

# Create new genotypes table
cur.execute("""CREATE TABLE isolates 
            (
              isolate TEXT,
              seg_1 TEXT,
              seg_2 TEXT,
              seg_3 TEXT,
              seg_4 TEXT,
              seg_5 TEXT,
              seg_6 TEXT,
              seg_7 TEXT,
              seg_8 TEXT,
              PRIMARY KEY(isolate)
            )""")

conn.commit()


cur.execute("SELECT isolate, primary_accession, segment from meta_data where isolate IS NOT NULL;")
rows = cur.fetchall()
# conn.close()

print("starting")
isolates_map = {}
for row in rows:
    isolate = row['isolate']
    if isolate not in isolates_map:
      isolates_map[isolate] = {1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None}
    if row['segment']:
      isolates_map[isolate][int(float(row['segment']))] = row['primary_accession']
print("inserting")
for key in isolates_map.keys():
  isolates = isolates_map[key]
  cur.execute("""
            INSERT INTO isolates (
              isolate,
              seg_1,
              seg_2,
              seg_3,
              seg_4,
              seg_5,
              seg_6,
              seg_7,
              seg_8
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            key,
            isolates[1],
            isolates[2],
            isolates[3],
            isolates[4],
            isolates[5],
            isolates[6],
            isolates[7],
            isolates[8]
        ))
  
conn.commit()
conn.close()


