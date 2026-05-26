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
from rest_framework import status
from Bio import SeqIO
from pathlib import Path



BASE_DIR = Path(__file__).resolve().parent.parent # goes up from api/
TASKS_DIR = BASE_DIR / "tasks"
RESOURCES_DIR = BASE_DIR / "resources"


def blastn(tmp_dir, query_path, results_path, db_path):
    
    if not os.path.exists(f"{db_path}.fa"):
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
                    queries[accession]["epa-ng"]["major_lwr"] = row['LWR']
                    break
        with open(minor_per_query_file, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['name'] == accession:
                    queries[accession]["epa-ng"]["minor"] = row['taxopath']
                    queries[accession]["epa-ng"]["minor_lwr"] = row['LWR']
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

    blast_db_path = BASE_DIR / "db" / "blast" / "db" 
    
    print("BLASTS_DB", blast_db_path)

    blastn(tmp_dir=tmp_dir, 
            query_path=query_path,
            results_path=results_path,
            db_path=blast_db_path)

    query_tophits_file = results_path / 'query_tophits.tsv'
    if not os.path.exists(query_tophits_file):
        raise ValueError("BLAST did not run.")
    
    blast_results = parse_query_tophits(query_tophits_file)

    # parse_blast_results(database=database, query_tophits=query_tophits_file, tmp_dir_input=inputs_path)

    ref_aln = RESOURCES_DIR / "reference_alignment.fa"
    ref_tree = RESOURCES_DIR / "taxon" / "ref_tree.treefile"
    taxon_major = RESOURCES_DIR / "taxon" / "taxon_major.tsv"
    taxon_minor = RESOURCES_DIR / "taxon" / "taxon_minor.tsv"
    meta_data = RESOURCES_DIR / "meta_data.tsv"


    aligned_out = str(results_path) + "/input_seqs_with_ref_alignment.fa"


    "ref_aln == reference sequences for all the reference sequences"
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

    results, tree = parse_clade_assignment_results(blast_results, results_path)



    return {"queries":results, "tree":tree}


