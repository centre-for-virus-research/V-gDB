import sqlite3

conn = sqlite3.connect("/Volumes/My Passport/CVR/gdb/RABV/rabv-gDB_Dec022025.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

# # Drop old table if it exists
cur.execute("DROP TABLE IF EXISTS genotypes;")

# Create new genotypes table
cur.execute("""CREATE TABLE genotypes 
            (
              major_clade TEXT,
              minor_clade TEXT
            )""")

conn.commit()


cur.execute("SELECT EPA_major_clade, EPA_minor_clade from meta_data WHERE EPA_major_clade IS NOT NULL;")
rows = cur.fetchall()
# conn.close()

print("starting")
major_clade_map = {}
for row in rows:
    major_clade = row['EPA_major_clade']
    minor_clade = row['EPA_minor_clade']
    if major_clade not in major_clade_map:
      major_clade_map[major_clade] = []
    if minor_clade:
      if minor_clade not in major_clade_map[major_clade]:
        major_clade_map[major_clade].append(minor_clade)
print("inserting")
print(major_clade_map)
for key in major_clade_map.keys():
  minor_clades = major_clade_map[key]
  if len(minor_clades) == 0:
    cur.execute("""
              INSERT INTO genotypes (
                major_clade, 
                minor_clade
              )
              VALUES (?, ?)
          """, (
              key,
              None,
          ))
  for minor_clade in minor_clades:
    cur.execute("""
              INSERT INTO genotypes (
                major_clade, 
                minor_clade
              )
              VALUES (?, ?)
          """, (
              key,
              minor_clade,
          ))
  
conn.commit()
conn.close()



# print("starting")
# isolates_map = {}
# for row in rows:
#     isolate = row['Parsed_Strain']
#     if isolate not in isolates_map:
#       isolates_map[isolate] = {1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None}
#     if row['segment']:
#       isolates_map[isolate][int(float(row['segment']))] = row['primary_accession']
# print("inserting")
# for key in isolates_map.keys():
#   isolates = isolates_map[key]
#   cur.execute("""
#             INSERT INTO isolates (
#               strain, 
#               seg_1,
#               seg_2,
#               seg_3,
#               seg_4,
#               seg_5,
#               seg_6,
#               seg_7,
#               seg_8
#             )
#             VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
#         """, (
#             key,
#             isolates[1],
#             isolates[2],
#             isolates[3],
#             isolates[4],
#             isolates[5],
#             isolates[6],
#             isolates[7],
#             isolates[8]
#         ))
  
# conn.commit()
# conn.close()

# 
