from django.db import connections
from os.path import join
import sqlite3
import os
import django
import sys
import sqlite3
from os.path import join
import subprocess
import csv
from models.clade_assignment.clade_assignment import CladeAssignment

from Bio import SeqIO
from pathlib import Path



# def extract_ref_seq():
#     print('Extracting Reference Sequences')
#     db_path = '/Users/danaallen/CVR/gdb/web-resources/V-gDB/db'
#     write_file = open(join(db_path, "db.fa"), 'w')
#     db_file = '/Volumes/My Passport/CVR/gdb/RABV/RABV-gDB_feb172026.db'
#     conn = sqlite3.connect(db_file)
#     cursor = conn.cursor()
#     cursor.execute("SELECT primary_accession, accession_type FROM meta_data where accession_type='reference' OR accession_type = 'master'")
#     ref_accs = cursor.fetchall()
#     ref_accs_list = [(item[0], item[1]) for item in ref_accs]
#     for each_acc, accession_type in ref_accs_list:
#         cursor.execute("SELECT sequence FROM sequences WHERE header = ?", [each_acc])
#         if each_acc == "NC_001542": accession_type="master"
#         result = cursor.fetchone()
#         if result:
#             sequence = result[0]
#             acc_type = "|" + accession_type
#             write_file.write(">" + each_acc + acc_type)
#             write_file.write("\n")
#             write_file.write(sequence)
#             write_file.write("\n")
#     write_file.close()
#     print('Finished creating blast database')

# def run_makeblastdb(db_fasta):
#     db_file_name = db_fasta
#     command = [
#         'makeblastdb',
#         '-in', db_fasta,
#         '-out', join(db_file_name),
#         '-title', "alignment",
#         '-dbtype', 'nucl'
#     ]
#     try:
#         subprocess.run(command, check=True)
#         print(f"makeblastdb ran successfully on {db_file_name}")
#     except subprocess.CalledProcessError as e:
#         print(f"Error running makeblastdb: {e}")






# def get_reference_alignment():

#     db_file = '/Volumes/My Passport/CVR/gdb/RABV/RABV-gDB_feb172026.db'
#     ref_out = '/Users/danaallen/CVR/gdb/web-resources/V-gDB/tmp/inputs/reference_alignment.fa'
   

#     conn = sqlite3.connect(db_file)
#     cursor = conn.cursor()

#     cursor.execute("""
#         SELECT sa.sequence_id, sa.alignment
#         FROM sequence_alignment sa
#         JOIN meta_data md ON md.primary_accession = sa.sequence_id
#         WHERE md.accession_type='reference'
#            OR md.accession_type='master'
#     """)

#     rows = cursor.fetchall()
#     # Get column names dynamically
#     col_names = [description[0] for description in cursor.description]

#      # Append to existing fasta file
#     with open(ref_out, 'a') as f:
#         for reference_accession,reference_alignment in rows:
#             f.write(f">{reference_accession}\n")
#             f.write(reference_alignment.strip() + "\n")

#     conn.close()



# def get_taxon_major_minor():

#     db_file = '/Volumes/My Passport/CVR/gdb/RABV/RABV-gDB_feb172026.db'
#     major_out = '/Users/danaallen/CVR/gdb/web-resources/V-gDB/tmp/inputs/taxon/taxon_major.tsv'
#     minor_out = '/Users/danaallen/CVR/gdb/web-resources/V-gDB/tmp/inputs/taxon/taxon_minor.tsv'

#     conn = sqlite3.connect(db_file)
#     cursor = conn.cursor()

#     cursor.execute("""
#         SELECT primary_accession, EPA_major_clade
#         FROM meta_data
#         WHERE (accession_type='reference' OR accession_type='master')
#         AND EPA_major_clade IS NOT NULL
#     """)


#     rows = cursor.fetchall()

#     # Write MAJOR file
#     with open(major_out, "w", newline="") as f_major:
#         writer = csv.writer(f_major, delimiter="\t")
#         writer.writerow(["acc", "sub_clade"])  # optional header
#         for acc, major in rows:
#             writer.writerow([acc, major])
        
#     cursor.execute("""
#         SELECT primary_accession, EPA_minor_clade
#         FROM meta_data
#         WHERE (accession_type='reference' OR accession_type='master')
#         AND EPA_minor_clade IS NOT NULL
#     """)


#     rows = cursor.fetchall()
#     # Write MINOR file
#     with open(minor_out, "w", newline="") as f_minor:
#         writer = csv.writer(f_minor, delimiter="\t")
#         writer.writerow(["acc", "sub_clade"])  # optional header
#         for acc, minor in rows:
#             writer.writerow([acc, minor])

#     conn.close()

# def get_reference_meta_data():

#     db_file = '/Volumes/My Passport/CVR/gdb/RABV/RABV-gDB_feb172026.db'
#     meta_out = '/Users/danaallen/CVR/gdb/web-resources/V-gDB/tmp/inputs/meta_data.tsv'
   

#     conn = sqlite3.connect(db_file)
#     cursor = conn.cursor()

#     cursor.execute("""
#         SELECT *
#         FROM meta_data
#         WHERE accession_type='reference'
#            OR accession_type='master'
#     """)

#     rows = cursor.fetchall()
#     # Get column names dynamically
#     col_names = [description[0] for description in cursor.description]

#     # Write everything to TSV
#     with open(meta_out, "w", newline="") as f:
#         writer = csv.writer(f, delimiter="\t")
#         writer.writerow(col_names)  # dynamic header
#         writer.writerows(rows)      # write all rows without unpacking

#     conn.close()


# def get_reference_tree():
#     db_file = '/Volumes/My Passport/CVR/gdb/RABV/RABV-gDB_feb172026.db'
#     tree_out = '/Users/danaallen/CVR/gdb/web-resources/V-gDB/tmp/inputs/taxon/ref_tree.treefile'

#     conn = sqlite3.connect(db_file)
#     cursor = conn.cursor()

#     cursor.execute("""
#                     SELECT newick
#                     FROM trees
#                     WHERE tree_type='reference'
#                 """)


#     row = cursor.fetchone()  # single row
#     if row:
#         newick_str = row[0]  # get the string from the tuple
#         with open(tree_out, "w") as f:
#             f.write(newick_str + "\n")  # write to .treefile
#     else:
#         print("No reference tree found in the database.")

#     conn.close()

BASE_DIR = Path(__file__).resolve().parent.parent # goes up from api/
TASKS_DIR = BASE_DIR / "tasks"
RESOURCES_DIR = BASE_DIR / "resources"


def blastn(tmp_dir, query_path, results_path, db_path):
    
    if not os.path.exists(db_path):
        raise ValueError("BLAST database not found, please message site administrator")
        
    db_file_name = os.path.basename(db_path)

    output_file = join(results_path,"query_tophits.tsv")
    command = [
        'blastn',
        '-query', query_path,
        '-db', join(db_path),
        '-task', 'blastn',
        '-max_target_seqs', '1',
        '-max_hsps', '1',
        '-out', output_file,
        '-outfmt', "6 qacc sacc pident sstrand"
    ]
    try:
        subprocess.run(command, check=True)
        print(f"blastn ran successfully. Results saved in {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running blastn: {e}")

def parse_blast_results(database, query_tophits, tmp_dir_input):
    with open(query_tophits, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            #query, ref, identity, strand,
            col1, col2, col3, col4 = row[0], row[1].split("|")[0], float(row[2]), row[3]
            combine_query_and_ref(database, col2, f'{tmp_dir_input}/{col1}.fa')

def combine_query_and_ref(database, reference_accession, query_path):

    # db_file = '/Volumes/My Passport/CVR/gdb/RABV/RABV-gDB_feb172026.db'

    # conn = sqlite3.connect(db_file)
    # cursor = conn.cursor()
    with connections[database].cursor() as cursor:
        cursor.execute(
            "SELECT alignment FROM sequence_alignment WHERE sequence_id = %s",
            [reference_accession]
        )

        result = cursor.fetchone()

    if result is None:
        raise ValueError(f"No alignment found for {reference_accession}")

    reference_alignment = result[0]

    # Append to existing fasta file
    with open(query_path, 'a') as f:
        f.write(f">{reference_accession}\n")
        f.write(reference_alignment.strip() + "\n")


def parse_query_tophits(query_tophits_file):

    results = {}
    with open(query_tophits_file, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            #query, ref, identity, strand,
            col1, col2, col3, col4 = row[0], row[1].split("|")[0], float(row[2]), row[3]
            results[col1] = {"blast_results":{
                                                "ref":col2,
                                                "identity":col3,
                                            }
                            }

    return results

def parse_clade_assignment_results(queries, results_dir):
    major_per_query_file = results_dir / "gappa_major_clades_assigned/per_query.tsv"
    minor_per_query_file = results_dir / "gappa_minor_clades_assigned/per_query.tsv"
    fasta_file = results_dir / "input_seqs_with_ref_alignment.fa"
    tree_file = results_dir / "gappa_major_clades_assigned/epa_result.newick"
    results = {}
    for accession in queries.keys():
        print(accession)
        queries[accession]["epa-ng"] = {}
        with open(major_per_query_file, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                
                if row['name'] == accession:
                    # print(row)
                    queries[accession]["epa-ng"]["major"] = row['taxopath']
                    break
        with open(minor_per_query_file, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['name'] == accession:
                    queries[accession]["epa-ng"]["minor"] = row['taxopath']
                    break
        # Iterate through the FASTA records
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id == accession:
                # print(f"ID: {record.id}")
                # print(f"Sequence: {record.seq}")
                queries[accession]["aligned_sequence"] = str(record.seq)
                break
        # else:
            # print("ID not found in FASTA.")

    # Read Newick as a single string
    with open(tree_file, "r") as f:
        tree = f.read().strip()

    return queries, tree 

def phylogenetic_clade_assignment_analysis(database, job_id):

    # Step 1: Run blast search: 
    tmp_dir = TASKS_DIR / job_id
    inputs_path = tmp_dir / "inputs" 
    query_path = inputs_path / "input.fa"
    results_path = tmp_dir / "results"
    

    blastn(tmp_dir=tmp_dir, 
            query_path=query_path,
            results_path=results_path,
            db_path='/Users/danaallen/CVR/gdb/web-resources/V-gDB/db/db.fa')

    query_tophits_file = results_path / 'query_tophits.tsv'
    if not os.path.exists(query_tophits_file):
        raise ValueError("BLAST did not run.")
    
    blast_results = parse_query_tophits(query_tophits_file)
    print(blast_results)

    # parse_blast_results(database=database, query_tophits=query_tophits_file, tmp_dir_input=inputs_path)

    ref_aln = RESOURCES_DIR / "reference_alignment.fa"
    ref_tree = RESOURCES_DIR / "taxon" / "ref_tree.treefile"
    taxon_major = RESOURCES_DIR / "taxon" / "taxon_major.tsv"
    taxon_minor = RESOURCES_DIR / "taxon" / "taxon_minor.tsv"
    meta_data = RESOURCES_DIR / "meta_data.tsv"


    aligned_out = str(results_path) + "/input_seqs_with_ref_alignment.fa"
    print(aligned_out)

    print("STARTING THE ANALYSIS")
    runner = CladeAssignment(ref_aln=ref_aln,
                            ref_tree=ref_tree,
                            query_fa=inputs_path/'input.fa',
                            taxon_major=taxon_major,
                            taxon_minor=taxon_minor,
                            meta_data=meta_data,
                            aligned_out=aligned_out,
                            output_dir=results_path
                            )

    runner.validate_inputs()

    try:
        runner.run_all()
    except subprocess.CalledProcessError as e:
        runner._die("Command failed with exit code " + str(e.returncode))


    # results = {"query":{},
    #         "tree":''}

    # 

    # print(results)


    results, tree = parse_clade_assignment_results(blast_results, results_path)



    return {"queries":results, "tree":tree}




def main():
    
    # extract_ref_seq()
    # run_makeblastdb("/Users/danaallen/CVR/gdb/web-resources/V-gDB/db/db.fa")
    # blastn(tmp_dir='/Users/danaallen/CVR/gdb/web-resources/V-gDB/tmp', 
    #         query_path='/Users/danaallen/CVR/gdb/web-resources/V-gDB/models/clade_assignment/query.fa',
    #         db_path='/Users/danaallen/CVR/gdb/web-resources/V-gDB/db/db.fa')
    # parse_blast_results(query_tophits='/Users/danaallen/CVR/gdb/web-resources/V-gDB/tmp/results/query_tophits.tsv')


    # get_taxon_major_minor()
    # get_reference_tree()
    # get_reference_meta_data()
    # get_reference_alignment()
    # ref_tree = get_reference_tree()
    # taxon_major = get_taxon_major()
    # taxon_minor = get_taxon_minor()




    run_phylogenetic_clade_assignment_analysis()

if __name__ == "__main__":
    main()