import sqlite3

conn = sqlite3.connect("/Volumes/My Passport/CVR/gdb/Flu/flu-gDB_dec02.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

# # Drop old table if it exists
cur.execute("DROP TABLE IF EXISTS trees;")

# Create new adaptive mutations table
cur.execute("""CREATE TABLE trees 
            (
              tree_name TEXT,
              tree_type TEXT,
              segment TEXT,
              newick TEXT
            )""")

conn.commit()
conn.close()
