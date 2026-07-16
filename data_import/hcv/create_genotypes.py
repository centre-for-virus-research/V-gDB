import sqlite3

conn = sqlite3.connect("/Volumes/My Passport/CVR/gdb/HCV/HCV_full.db")
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


cur.execute("SELECT nearest_reference_genotype as EPA_major_clade, nearest_reference_subtype as EPA_minor_clade from meta_data WHERE nearest_reference_genotype IS NOT NULL;")
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
    
cur.execute(f"""
            ALTER TABLE meta_data
            RENAME COLUMN nearest_reference_subtype TO EPA_minor_clade;"""
            )
cur.execute(f"""
            ALTER TABLE meta_data
            RENAME COLUMN nearest_reference_genotype TO EPA_major_clade;"""
            )
  
conn.commit()
conn.close()


