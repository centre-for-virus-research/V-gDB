from argparse import ArgumentParser
import sqlite3
import csv

# "/Volumes/My Passport/V-gDB-flu-Sep292025.db"
# "/Volumes/My Passport/CVR/gdb/Flu/adaptive_mutations/Clean PB2.csv"
def import_adaptive_mutations(file, database):
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row  
    cur = conn.cursor()

    with open(file, newline="") as csvfile:
        print(f"Starting import for {file}")

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

        print(f"Finished adaptive mutations import for {file}")
    
    conn.commit()
    conn.close()



if __name__ == "__main__":
    parser = ArgumentParser(description='Imports adaptive mutations from CSV file')
    parser.add_argument('-f', '--file', help='adaptive mutation (.csv) file', required=True)
    parser.add_argument('-d', '--database', help='path to database (.db) file', required=True)
    args = parser.parse_args()
    
    import_adaptive_mutations(args.file, args.database)
