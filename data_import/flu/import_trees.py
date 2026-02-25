from argparse import ArgumentParser
import sqlite3
import csv

# /Volumes/My Passport/CVR/gdb/Flu/flu-gDB_dec02.db
"/Volumes/My Passport/CVR/gdb/Flu/adaptive_mutations/Clean PB2.csv"

"/Volumes/My Passport/CVR/gdb/Flu/trees/seg1_cluster_rep.nwk"
def import_trees(file, database):
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row  
    cur = conn.cursor()

    with open("/Volumes/My Passport/CVR/gdb/Flu/trees/seg1_cluster_rep.nwk", "r") as f:
        data = f.read()
        cur.execute("""
                INSERT INTO trees (
                tree_name,
                tree_type,
                segment,
                newick
                )
                VALUES (?, ?, ?, ?)
            """, (
                "Segment 1 - Cluster Representative",
                "cluster",
                "1",
                data
            ))
            

        print(f"Finished adaptive mutations import for {file}")
    
    conn.commit()
    conn.close()



if __name__ == "__main__":
    parser = ArgumentParser(description='Imports trees')
    parser.add_argument('-f', '--file', help='newick tree (.nwk) file', required=True)
    parser.add_argument('-d', '--database', help='path to database (.db) file', required=True)
    args = parser.parse_args()
    
    import_trees(args.file, args.database)
