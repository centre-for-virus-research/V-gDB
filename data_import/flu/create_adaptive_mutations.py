import sqlite3

conn = sqlite3.connect("/Volumes/My Passport/V-gDB-flu-Sep292025.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

# # Drop old table if it exists
cur.execute("DROP TABLE IF EXISTS adaptive_mutations;")

# Create new adaptive mutations table
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
conn.close()
