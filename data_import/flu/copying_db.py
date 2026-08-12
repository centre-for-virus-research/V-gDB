import sqlite3

# Path to the source and destination databases
source_db = "/Volumes/My Passport/V-gDB-flu-Sep292025.db"
dest_db = "/Volumes/My Passport/CVR/gdb/Flu/flu-gDB_dec02.db"
table_name = "features"

# Connect to the destination database
conn = sqlite3.connect(dest_db)
cur = conn.cursor()



# Attach the source database
cur.execute(f"ATTACH DATABASE '{source_db}' AS src_db")

# Create the table in the destination if it doesn't exist
# (You can skip this if the table structure already exists)
cur.execute(f"""
    CREATE TABLE IF NOT EXISTS {table_name} AS
    SELECT * FROM src_db.{table_name} WHERE 0
""")

# Copy all data from source to destination
cur.execute(f"""
    INSERT INTO {table_name}
    SELECT * FROM src_db.{table_name}
""")

# Commit and close
conn.commit()
conn.close()

print(f"Table '{table_name}' copied from {source_db} to {dest_db}")