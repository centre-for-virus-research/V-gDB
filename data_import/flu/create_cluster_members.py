import sqlite3

conn = sqlite3.connect("/Volumes/My Passport/CVR/gdb/Flu/flu-gDB_dec02.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

# # Drop old table if it exists
cur.execute("DROP TABLE IF EXISTS cluster_members;")

# Create new adaptive mutations table
cur.execute("""CREATE TABLE cluster_members 
            (
              segment TEXT,
              primary_accession TEXT
            )""")

conn.commit()
conn.close()

print("cluster_members was successfully created")
