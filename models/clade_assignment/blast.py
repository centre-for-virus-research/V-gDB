#python blastAlignment.py -s Y -f generic-influenza/ref_list.txt

import os
import csv
import sys
import shutil
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
from ete3 import Tree

def prune_tree():

    tree = Tree("/Users/danaallen/CVR/gdb/web-resources/V-gDB/resources/HCV/taxon/ref_tree.treefile")

    msa_taxa = {record.id for record in SeqIO.parse("/Users/danaallen/CVR/gdb/web-resources/V-gDB/resources/HCV/reference_alignment.fa", "fasta")}
    tree_names = {leaf.name for leaf in tree}   

    common = msa_taxa & tree_names


    tree.prune(common, preserve_branch_length=True)
    # tree.write(outfile="pruned.nwk")

    # tree.prune(msa_taxa, preserve_branch_length=True)

    tree.write(outfile="/Users/danaallen/CVR/gdb/web-resources/V-gDB/resources/HCV/taxon/ref_tree.treefile")

def add_reference_header():
    records = []
    for record in SeqIO.parse("/Users/danaallen/CVR/gdb/web-resources/V-gDB/db/HCV/db.fa", "fasta"):
        record.id = f"{record.id}|reference"
        record.description = record.id
        records.append(record)

    SeqIO.write(records, "/Users/danaallen/CVR/gdb/web-resources/V-gDB/db/HCV/db.fa", "fasta")


def run_makeblastdb():
    tmp_dir = "/Users/danaallen/CVR/gdb/web-resources/V-gDB/db"
    os.makedirs(join(tmp_dir, 'HCV'), exist_ok=True)
    db_fasta = "/Users/danaallen/CVR/gdb/web-resources/V-gDB/db/HCV/db.fa"

    command = [
        'makeblastdb',
        '-in', db_fasta,
        '-out', join(tmp_dir, 'HCV', "db"),
        '-title', "alignment",
        '-dbtype', 'nucl'
    ]
    try:
        subprocess.run(command, check=True)
        print(f"makeblastdb ran successfully on {db_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running makeblastdb: {e}")

def process():
    # add_reference_header()
    # run_makeblastdb()
    prune_tree()
				

if __name__ == "__main__":
	
	process()