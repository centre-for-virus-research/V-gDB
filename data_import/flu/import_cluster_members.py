from argparse import ArgumentParser
import sqlite3
import csv

# "/Volumes/My Passport/V-gDB-flu-Sep292025.db"
"/Volumes/My Passport/CVR/gdb/Flu/adaptive_mutations/Clean PB2.csv"
def import_cluster_members(file, database, segment):
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row  
    cur = conn.cursor()

    unique_clusters = set()
    with open(file, newline="") as csvfile:
        print(f"Starting import for {file}")

        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if not row:          # skip empty rows
                continue
            unique_clusters.add(row[0])
        cur.executemany(
            """
            INSERT INTO cluster_members (segment, primary_accession)
            VALUES (?, ?)
            """,
            [(segment, value) for value in unique_clusters]
        )

        print(f"Finished cluster_member import for {file}")
    
    conn.commit()
    conn.close()



if __name__ == "__main__":
    parser = ArgumentParser(description='Imports cluster_members from CSV file')
    parser.add_argument('-f', '--file', help='cluster segment (.csv) file', required=True)
    parser.add_argument('-d', '--database', help='path to database (.db) file', required=True)
    parser.add_argument('-s', '--segment', help='segment', required=True)
    args = parser.parse_args()
    
    import_cluster_members(args.file, args.database, args.segment)
