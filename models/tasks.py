from django.db import connections
import time
import logging
import django_rq
import django
from itertools import groupby
import rq

import os
import sys
import re
import csv
import uuid
import glob
import shutil
import sqlite3
import requests
import subprocess
from Bio import SeqIO
from time import sleep
from Bio.Seq import Seq
from os.path import join
from argparse import ArgumentParser
from pathlib import Path
# from FeatureCalculator  import FeatureCordCalculator 


project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'backend'))
sys.path.insert(0, project_root)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'backend.settings')

# Initialize Django
django.setup()
def fasta(fasta_name):
	read_file = open(fasta_name)
	faiter = (x[1] for x in groupby(read_file, lambda line: line[0] == ">"))
	for header in faiter:
		#header = header.next()[1:].strip()
		header = next(header)[1:].strip()
		seq = "".join(s.strip() for s in next(faiter))
		yield header, seq
		
sequences ='>MT895968' \
'acgcttaacaaccagatcaaagaaaaaacagacattgtcaattgcaaagcaaaactgtaacacccctacaatggatgccgacaagattgtattcaaagtcaataatcaggtggtctctttgaagcctgagattatcgtggatcaacatgagtacaagtaccctgccatcaaagatttgaaaaagccctgtataaccctaggaaaggctcccgatttaaataaagcatacaagtcagttttgtcaggcatgagcgccgccaaacttgatcctgacgatgtatgttcctatttggcagcggcaatgcagttttttgaggggacatgtccggaagactggaccagctatggaatcgtgattgcacgaaaaggagataagatcaccccaggttctctggtggagataaaacgtactgatgtagaagggaattgggctctgacaggaggcatggaactgacaagagaccccactgtccctgagcatgcgtccttagtcggtcttctcttgagtctgtataggttgagcaaaatatccgggcaaaacactggtaactataagacaaacattgcagacaggatagagcagatttttgagacagccccttttgttaaaatcgtggaacaccatactctaatgacaactcacaaaatgtgtgctaattggagtactataccaaacttcagatttttggccggaacctatgacatgtttttctcccggattgagcatctatattcagcaatcagagtgggcacagttgtcactgcttatgaagactgttcaggactggtatcatttactgggttcataaaacaaatcaatctcaccgctagagaggcaatactatatttcttccacaagaactttgaggaagagataagaagaatgtttgagccagggcaggagacagctgttcctcactcttatttcatccacttccgttcactaggcttgagtgggaaatctccttattcatcaaatgctgttggtcacgtgttcaatctcattcactttgtaggatgctatatgggtcaagtcagatccctaaatgcaacggttattgctgcatgtgctcctcatgaaatgtctgttctagggggctatctgggagaggaattcttcgggaaagggacatttgaaagaagattcttcagagatgagaaagaacttcaagaatacgaggcggctgaactgacaaagactgacgtagcactggcagatgatggaactgtcaactctgacgacgaggactactgctcaggtgaaaccagaagtccggaggctgtttatactcgaatcatgatgaatggaggtcgactaaagagatctcacatacggagatatgtctcagtcagttccaatcatcaagcccgtccaaactcattcgccgagtttctaaacaagacatattcgagtgactcataagaagttgaataacaaaatgccggaaatctacggattgtgtatatccatcatgaaaaaaactaacacccctcctttcgaaccatcccaaacatgagcaagatctttgtcaatcctagtgctattagagccggtctggccgatcttgagatggctgaagaaactgttgatctgatcaatagaaatatcgaagacaatcaggctcatctccaaggggaacccatagaagtggacaatctccctgaggatatggggcgacttcacctggatgatggaaaatcgcccaaccctggtgagatggccaaggtgggagaaggcaagtatcgagaggactttcagatggatgaaggagaggatcctagcttcctgttccagtcatacctggaaaatgttggagtccaaatagtcagacaaatgaggtcaggagagagatttctcaagatatggtcacagaccgtagaagagattatatcctatgtcgcggtcaactttcccaaccctccaggaaagtcttcagaggataaatcaacccagactactggccgagagctcaagaaggagacaacacccactccttctcagagagaaagccaatcatcgaaagccaggatggcggctcaaactgcttctggccctccagcccttgaatggtcggccgccaatgaagaggatgatctatcagtggaggctgagatcgctcaccagattgcagaaagtttctccaaaaaatataagtttccctctcgatcctcagggatactcttgtataattttgagcaattgaaaatgaaccttgatgatatagttaaagaggcaaaaaatgtaccaggtgtgacccgtttagcccatgacgggtccaaactccccctaagatgtgtactgggatgggtcgctttggccaactctaagaaattccagttgttagtcgaatccgacaagctgagtaaaatcatgcaagatgacttgaatcgctatacatcttgctaaccgagcctctccactcagtccctctagacaataaagtccgagatgtcctaaagtcaacatgaaaaaaacaggcaacaccactgataaaatgaactttctacgtaagatagtgaaaaattgcagggacgaggacactcaaaaaccctctcccgtgtcagcccctctggatgacgatgacttgtggcttccaccccctgaatacgtcccgctgaaagaacttacaagcaagaagaacatgaggaacttttgtatcaacggaggggttaaagtgtgtagcccgaatggttactcgttcaggatcctgcggcacattctgaaatcattcgacgagatatattctgggaatcataggatgatcgggttagtcaaagtagttattggactggctttgtcaggatctccagtccctgagggcatgaactgggtatacaaattgaggagaacctttatctttcagtgggctgattccaggggccctcttgaaggggaggagttggaatactctcaggagatcacttgggatgatgatactgagttcgtcggattgcaaataagagtgattgcaaaacagtgtcatatccagggcagaatctggtgtatcaacatgaacccgagagcatgtcaactatggtctgacatgtctcttcagacacaaaggtccgaagaggacaaagattcctctctgcttctagaataatcagattatatcccgcaaatttatcacttgtttacctctggaggagagaacatatgggctcaactccaacccttgggagcaatataacaaaaaacatgttatggtgccattaaaccgctgcatttcatcaaagtcaagttgattacctttacattttgatcctcttggatgtgaaaaaaactattaacatccctcaaaagactcaaggaaagatggttcctcaggctctcctgtttgtaccccttctggtttttccattgtgttttgggaaattccctatttacacgataccagacaagcttggtccctggagcccgattgacatacatcacctcagctgcccaaacaatttggtagtggaggacgaaggatgcaccaacctgtcagggttctcctacatggaacttaaagttggatacatcttagccataaaaatgaacgggttcacttgcacaggcgttgtgacggaggctgaaacctacactaacttcgttggttatgtcacaaccacgttcaaaagaaagcatttccgcccaacaccagatgcatgtagagccgcgtacaactggaagatggccggtgaccccagatatgaagagtctctacacaatccgtaccctgactaccgctggcttcgaactgtaaaaaccaccaaggagtctctcgttatcatatctccaagtgtggcagatttggacccatatgacagatcccttcactcgagggtcttccctagcgggaagtgctcaggagtagcggtgtcttctacctactgctccactaaccacgattacaccatttggatgcccgagaatccgagactagggatgtcttgtgacatttttaccaatagtagagggaagagagcatccaaagggagtgagacttgcggctttgtagatgaaagaggcctatataagtctttaaaaggagcatgcaaactcaagttatgtggagttctaggacttagacttatggatggaacatgggtcgcgatgcaaacatcaaatgaaaccaaatggtgccctcccgatcagttggtgaacctgcacgactttcgctcagacgaaattgagcaccttgttgtagaggagttggtcaggaagagagaggagtgtctggatgcactagagtccatcatgacaaccaagtcagtgagtttcagacgtctcagtcatttaagaaaacttgtccctgggtttggaaaagcatataccatattcaacaagaccttgatggaagccgatgctcactacaagtcagtcagaacttggaatgagatcctcccttcaaaagggtgtttaagagttggggggaggtgtcatcctcatgtgaacggggtgtttttcaatggtataatattaggacctgacggcaatgtcttaatcccagagatgcaatcatccctcctccagcaacatatggagttgttggaatcctcggttatcccccttgtgcaccccctggcagacccgtctaccgttttcaaggacggtgacgaggctgaggattttgttgaagttcaccttcccgatgtgcacaatcaggtctcaggagttgacttgggtctcccgaactgggggaagtatgtattactgagtgcaggggccctgactgccttgatgttgataattttcctgatgacatgttgtagaagagtcaatcgatcagaacctacgcaacacaatctcagagggacagggagggaggtgtcagtcactccccaaagcgggaagatcatatcttcatgggaatcacacaagagtgggggtgagaccagactgtgaggactggccgtcctttcaacgatccaagtcctgaagatcacctccccttggggggttctttttgaaaaaaaacctgggttcaatagtcctccttgaactccatgcaactgggtagattcaagagtcatgagattttcattaatcctctcagttgatcaagcaagatcatgtagattctcataataggggagatcttctagcagtttcagtgactaacggtactttcattctccaggaactgacaccaacagttgtagacaaaccacggggtgtctcgggtgactctgtgcttgggcacagacaaaggtcatggtgtgttccatgatagcggactcaggatgagttaattgagagaggcagtcttcctcccgtgaaggacataagcagtagctcacaatcatctcgcgtctcagcaaagtgtgcataattataaagtgctgggtcatctaagcttttcagtcgagaaaaaaacattagatcagaagaacaactggcaacacttctcaacctgagacctacttcaagatgctcgatcctggagaggtctatgatgaccctattgacccaatcgagttagaggctgaacccagaggaacccccactgtccccaacatcttgaggaactctgactacaatctcaactctcctttgatagaagatcctgctagactaatgttagaatggttaaaaacagggaatagaccttatcggatgactctaacagacaattgctccaggtctttcagagttttgaaagattatttcaagaaggtagatttgggttctctcaaggtgggcggaatggctgcacagtcaatgatttctctctggttatatggtgcccactctgaatccaacaggagccggagatgtataacagacttggcccatttctattccaagtcgtcccccatagagaagctgttgaatctcacgctaggaaatagagggctgagaatccccccagagggagtgttaagttgccttgagagggttgattatgataatgcatttggaaggtatcttgccaacacgtattcctcttacttgttcttccatgtaatcaccttatacatgaacgccctagactgggatgaagaaaagaccatcctagcattatggaaagatttaacctcagtggacatcgggaaggacttggtaaagttcaaagaccaaatatggggactgctgatcgtgacaaaggactttgtttactcccaaagttccaattgtctttttgacagaaactacacacttatgctaaaagatcttttcttgtctcgcttcaactccttaatggtcttgctctctcccccagagccccgatactcagatgacttgatatctcaactatgccagctgtacattgctggggatcaagtcttgtctatgtgtggaaactccggctatgaagtcatcaaaatattggagccatatgtcgtgaatagtttagtccagagagcagaaaagtttaggcctctcattcattccttgggagactttcctgtatttataaaagacaaggtaagtcaacttgaagagacgttcggtccctgtgcaagaaggttctttagggctctggatcaattcgacaacatacatgacttggtttttgtgtatggctgttacaggcattgggggcacccatatatagattatcgaaagggtctgtcaaaactatatgatcaggttcacattaaaaaagtgatagataagtcctaccaggagtgcttagcaagcgacctagccaggaggatccttagatggggttttgataagtactccaagtggtatctggattcaagattcctagcccgagaccaccccttgactccttatatcaaaacccaaacatggccacccaaacatattgtagacttggtgggggatacatggcacaagctcccgatcacgcagatctttgagattcctgaatcaatggatccgtcagaaatattggatgacaaatcacattctttcaccagaacgagactagcttcttggctgtcagaaaaccgaggggggcctgttcctagcgaaaaagttattatcacggccctgtctaagccgcctgtcaatccccgagagtttctgaggtctatagacctcggaggattgccagatgaagacttgataattggcctcaagccaaaggaacgggaattgaagattgaaggtcgattctttgctctaatgtcatggaatctaagattgtattttgtcatcactgaaaaactcttggccaactacatcttgccactttttgacgcgctgactatgacagacaacctgaacaaggtgtttaaaaagctgatcgacagggtcaccgggcaagggcttttggactattcaagggtcacatatgcatttcacctggactatgaaaagtggaacaaccatcaaagattagagtcaacagaggatgtattttctgtcctagatcaagtgtttggattgaagagagtgttttctagaacacacgagttttttcaaaaggcctggatctattattcagacagatcagacctcatcgggttacgggaggatcaaatatactgcttagatgcgtccaacggcccaacctgttggaatggccaggatggcgggctagaaggcttacggcagaagggctggagtctagtcagcttattgatgatagatagagaatctcaaatcaggaacacaagaaccaaaatactagctcaaggagacaaccaggttttatgtccgacatatatgttgtcgccagggctatctcaagaggggctcctctatgaattggagagaatatcaaggaatgcactttcgatatacagagccgtcgaggaaggggcatctaagctagggctgatcatcaagaaagaagagaccatgtgtagttatgacttcctcatctatggaaaaacccctttgtttagaggtaacatattggtgcctgagtccaaaagatgggccagagtctcttgcgtctctaatgaccaaatagtcaacctcgccaatataatgtcgacagtgtccaccaatgcgctaacagtggcacaacactctcaatctttgatcaaaccgatgagggattttctgctcatgtcagtacaggcagtctttcactacctgctatttagcccaatcttaaagggaagagtttacaagattctgagcgctgaaggggagagctttctcctagccatgtcaaggataatctatctagatccttctttgggaggggtatctggaatgtccctcggaagattccatatacgacagttctcagaccctgtctctgaagggttatccttctggagagagatctggttaagctcccacgagtcctggattcacgcgttgtgtcaagaggctggaaacccagatcttggagagagaacactcgagagcttcactcgccttctagaagatcctaccaccttaaatatcagaggaggggccagtcctaccattctactcaaggatgcaatcagaaaggctttatatgacgaggtggacaaggtggagaattcagagtttcgagaggcaatcctgttgtccaagacccatagagataattttatactcttcttaacatctgttgagcctctgtttcctcgatttctcagtgagctattcagttcgtcttttttgggaatccccgagtcaatcattggattgatacaaaactcccgaacgataagaaggcagtttaaaaagagtctctcaaaaactttagaagagtccttctacaactcagagatccacgggattagtcggatgacccagacacctcagagggttgggggggtgtggccttgctcttcagagagggcagatctacttagggagatctcttggggaagaaaagtggtaggcacgacagttcctcacccttctgagatgttggggttacttcccaagtcctctatttcttgcacttgtggagcaacaggaggaggcaatcctagagtttctgtatcagtactcccgtcctttgatcagtcatttttttcacgaggccccctaaaggggtacttgggctcgtccacctctatgtcgacccagctattccatgcatgggaaaaagtcactaatgttcatgtggtgaagagagctctatcgttaaaagaatctataaactggttcattactagagattccaacttggctcaagctctaattaggaacattatgtctctgacaggccctgatttccctctagaggaggcccctgtcttcaaaaggacggggtcagccttgcataggttcaagtctgccagatacagcgaaggagggtattcttctgtctgcccgaacctcctctctcatatttctgttagtacagacaccatgtctgatttgacccaagacgggaagaactacgatttcatgttccagccattgatgctttatgcacagacatggacatcagagctggtacagagagacacaaggctaagagactctacgtttcattggcacctccgatgcaacaggtgtgtgagacccattgacgacgtgaccctggagacctctcagatcttcgagtttccggatgtgtcgaaaagaatatccagaatggtttctggggctgtgcctcacttccagaggcttcccgatatccgtctgagaccaggagattttgaatctctaagcggtagagaaaagtctcaccatatcggatcagctcaggggctcttatactcaatcttagtggcaattcacgactcaggatacaatgatggaaccatcttccctgtcaacatatacggcaaggtttcccctagagactatttgagagggctcgcaaggggagtattgataggatcctcgatttgcttcttgacaagaatgacaaatatcaatattaatagacctcttgaattgatctcaggggtaatctcatatattctcctgaggctagataaccatccctccttgtacataatgctcagagaaccgtctcttagaggagagatattttctatccctcagaaaatccccgccgcttatccaaccactatgaaagaaggcaacagatcaatcttgtgttatctccaacatgtgctacgctatgagcgagagataatcacggcgtctccagagaatgactggctatggatcttttcagactttagaagtgccaaaatgacgtacctaaccctcattacttaccagtctcatcttctactccagagggttgagagaaacatatctaagagtatgagagataacctgcgacaattgagttccttgatgaggcaggtgctgggcgggcacggagaagataccttagagtcagacgacaacattcaacgactgctaaaagactctttacgaaggacaagatgggtggatcaagaggtgcgccatgcagctagaaccatgactggagattacagccccaacaagaaggtgtcccgtaaggtaggatgttcagaatgggtctgctctgctcaacaggttgcagtctctacctcagcaaacccggcccctgtctcggagcttgacataagggccctctctaagaggttccagaaccctttgatctcgggcttgagagtggttcagtgggcaaccggtgctcattataagcttaagcctattctagatgatctcaatgttttcccatctctctgccttgtagttggggacgggtcaggggggatatcaagggcagtcctcaacatgtttccagatgccaagcttgtgttcaacagtcttttagaggtgaatgacctgatggcttccggaacacatccactgcctccttcagcaatcatgaggggaggaaatgatatcgtctccagagtgatagattttgactcaatctgggaaaaaccgtccgacttgagaaacttggcaacctggaaatacttccagtcagtccaaaagcaggtcaacatgtcctatgacctcattatttgcgatgcagaagttactgacattgcatctatcaaccggataaccctgttaatgtccgattttgcattgtctatagatggaccactctatttggtcttcaaaacttatgggactatgctagtaaatccaaactacaaggctattcaacacctgtcaagagcgttcccctcggtcacagggtttatcacccaagtaacttcgtctttttcatctgagctctacctccgattctccaaacgagggaagtttttcagagatgctgagtacttgacctcttccacccttcgagaaatgagccttgtgttattcaattgtagcagccccaagagtgagatgcagagagctcgttccttgaactgtcaggatcttgtgagaggatttcctgaagaaatcatatcaaatccttacaatgagatgatcataactntgattgacagtgatgtagaatcttttctagtccacaagatggttgatgatcttgagttacagaggggaactctgtctaaagtggctatcattatagccatcatgatagttttctccaacagagtcttcaacgtttccaaacccctaactgaccccttgttctatccaccgtctgatcccaaaatcctgaggcacttcaacatatgttgcagtactatgatgtatctatctactgctttaggtgacgtccctagcttcgcaagacttcacgacctgtataacagacctataacttattacttcagaaagcaagtcattcgagggaacgtttatctatcttggagttggtccaacgacacctcagtgttcaaaagggtagcctgtaattctagcctgagtctgtcatctcactggatcaggttgatttacaagatagtgaagactaccagactcgttggcagcatcaaggatctatccagagaagtggaaagacaccttcataggtacaacaggtggatcaccctagaggatatcagatctagatcatccctactagactacagttgcctgtgaaccggatactcctggaagcctgcccatgctaagactcttgtgtgatgtatcttgaaaaaaacaagatcctaaatctgaacctttggttgtttgattgtttttctcatttttgttgtttatttgttaagcgt'


def write_fasta_file(sequence_data, file_path):
    """
    Writes a sequence string in FASTA format to a file.

    Args:
    - sequence_data (str): The input sequence in ">id\nsequence" format.
    - file_path (str): The path where the FASTA file should be saved.
    """

    # lines = sequence_data.strip().split("\n")  # Ensure no extra spaces or newlines
    # if not lines or not lines[0].startswith(">"):
    #     raise ValueError("Invalid FASTA format: Missing sequence ID line")

    # fasta_id = lines[0]  # Extract the sequence ID (with ">")
    # sequence = "".join(lines[1:])  # Join sequence lines in case they are split

    # # Ensure line wrapping at 60 characters per line
    # formatted_sequence = "\n".join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
    fasta_sequence = sequence_data.replace("%0A", "\n").replace("%20", " ").replace("%3E", ">")

    # Write to file
    with open(file_path, "w") as fasta_file:
        # fasta_file.write(f"{fasta_id}\n{formatted_sequence}\n")
        fasta_file.write(fasta_sequence + "\n")

    print(f"FASTA file saved to {file_path}")
logger = logging.getLogger(__name__)
class Tasks:
    def __init__(self, database):
        self.database = database  

    @staticmethod
    def path_to_basename(file_path):
        path = os.path.basename(file_path)
        return path.split('.')[0]

    def example_task():
        job = rq.get_current_job()
        """A simple task that sleeps for n seconds."""
        logger.info(f"Starting the tasks")
        
        # time.sleep(5)
        logger.info(f"after first")
        job.meta['status'] = "running"
        
        logger.info(f"After processing")
        job.save_meta()
        time.sleep(3)
        job.meta['status'] = "running"
        logger.info(f"updating")
        job.save_meta()
        time.sleep(3)
        job.meta['status'] = "failed"
        logger.info(f"finished")
        job.save_meta()
        time.sleep(3)
        logger.info(f"done")
        job.meta['status'] = "done"
        job.save_meta()
        return "WOOOO WE ARE DONE!"
    
    def process_genbank_submission(self):
        db_path = Path(join("analysis","db.fa"))
        job = "test"
        os.makedirs(join("jobs",job), exist_ok=True)

    #     # write_fasta_file(sequences, "jobs/"+job+"/input.fasta")


        tmp_dir = "jobs/"+job
        query_path = join(tmp_dir, "input.fasta")


        os.makedirs(join(tmp_dir, 'analysis'), exist_ok=True)
        os.makedirs(join(tmp_dir, 'sorted_fasta'), exist_ok=True)
        os.makedirs(join(tmp_dir, 'merged_fasta'), exist_ok=True)
        os.makedirs(join(tmp_dir, 'sorted_all'), exist_ok=True)
        grouped_fasta = join(tmp_dir, 'grouped_fasta')
        os.makedirs(grouped_fasta, exist_ok=True)
        os.makedirs(join(tmp_dir, 'reference_sequences'), exist_ok=True)
        query_aln_output_dir = join(tmp_dir, "query_aln")
        os.makedirs(join(tmp_dir, 'master_sequences'), exist_ok=True)
        master_and_reference_merged = join(tmp_dir,  "master_and_reference_merged")
        mafft_reference_alignment = join(tmp_dir, "mafft_reference_alignment")
        reference_alignments =  join(tmp_dir, "reference_alignments")
        query_ref_alignment = join(tmp_dir, "query_ref_alignment")
        table2asn_tmp = join(tmp_dir, "Table2asn", "tmp")
        create_tmp_dir = os.makedirs(table2asn_tmp, exist_ok=True)

        os.makedirs(master_and_reference_merged, exist_ok=True)
        os.makedirs(mafft_reference_alignment, exist_ok=True)
        os.makedirs(reference_alignments, exist_ok=True)
        os.makedirs(query_ref_alignment, exist_ok=True)
                    
        if not db_path.is_file():
            self.__extract_ref_seq(db_path)

        self.blast_analysis(tmp_dir) 
        for each_ref_merged_acc in os.listdir(master_and_reference_merged):
            self.mafft_ref_sequences(join(master_and_reference_merged, each_ref_merged_acc), mafft_reference_alignment)

        for each_ref_aln in os.listdir(mafft_reference_alignment):
            output_file = self.path_to_basename(each_ref_aln)
            output_file = join(reference_alignments, output_file + '.fasta')
            self.extract_matching_sequences(join(mafft_reference_alignment), output_file)

        for each_query_seq in os.listdir(grouped_fasta):
            reference_file = self.path_to_basename(each_query_seq)
            self.mafft_query_sequences(join(grouped_fasta, each_query_seq), join(reference_alignments, reference_file), query_ref_alignment)

        
        return
    

    def blast_analysis(self, tmp_dir):

        values = {}
        records = {}

        query_path = join(tmp_dir, "input.fasta")
        ref_path = join(tmp_dir, "ref.fa")
        db_path = join("analysis", "db.fa")
        query_tophits = join(tmp_dir, "query_tophits.tsv")
        query_tophits_uniq = join(tmp_dir, "query_tophits_uniq.tsv")
        sorted_fasta = join(tmp_dir, "sorted_fasta")
        merged_fasta = join(tmp_dir, "merged_fasta")
        sorted_all = join(tmp_dir, "sorted_all")
        grouped_fasta = join(tmp_dir, "grouped_fasta")
        ref_seq_dir = join(tmp_dir, "reference_sequences")
        master_seq_dir = join(tmp_dir, 'master_sequences')
        master_and_reference_merged = join(tmp_dir, 'master_and_reference_merged')

        self.__run_makeblastdb(db_path)
        self.blastn(tmp_dir, query_path, db_path)

        with open(query_tophits, newline='') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                #query, ref, identity, strand,
                col1, col2, col3, col4 = row[0], row[1].split("|")[0], float(row[2]), row[3]
                if col1 in values:
                    existing_value = values[col1]
                    if col3 > existing_value:
                        records[col1] = [col1, col2, col3, col4]
                        values[col1] = col3
                else:
                    records[col1] = [col1, col2, col3, col4]
                    values[col1] = col3

        with open(query_tophits_uniq, 'w', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            for record in records.values():
                writer.writerow(record)

        seq_dicts = {}
        query_seqs = fasta(query_path)
        for rows in query_seqs:
            seq_dicts[rows[0].strip()] = rows[1].strip()

        #Seperate plus and minus strand sequences
        for each_line in open(query_tophits_uniq, 'r'):
            query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
            if strand == "plus":
                with open(join(sorted_fasta, "plus.fa"), 'a') as file_plus:
                    file_plus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")
            else:
                with open(join(sorted_fasta, "minus.fa"), 'a') as file_minus:
                    file_minus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")

        for each_file in os.listdir(sorted_fasta):
            if "minus" in each_file:
                command = ["seqkit", "seq", "-r", "-p", "-v", "-t", "dna", join(sorted_fasta, each_file), ">", join(merged_fasta, each_file)]
            else:
                command = ["cp", join(sorted_fasta, each_file), join(merged_fasta, each_file)]

            try:
                print(' '.join(command))
                os.system(' '.join(command))
                print(f"seqkit ran successfully for {each_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error running seqkit: {e}")

        file_list = []
        for each_file in os.listdir(merged_fasta):
            prefix = merged_fasta + "/"
            file_list.append(prefix + each_file)

        command = ["cat", " ".join(file_list), ">", join(sorted_all, "query_seq.fa")]
        print(' '.join(command))
        try:
            os.system(' '.join(command))
            print(f"concatenation sucessful")
        except subprocess.CalledProcessError as e:
            print(f"Error in concatenation: {e}")

        grouped_dict = {}
        for each_line in open(query_tophits_uniq, 'r'):
            query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
            if ref_acc not in grouped_dict:
                grouped_dict[ref_acc] = [query_acc]
            else:
                grouped_dict[ref_acc].append(query_acc)

        seq_dicts = {}
        master_seq_dict = {}
        query_seqs = fasta(join(sorted_all, "query_seq.fa"))
        for rows in query_seqs:
            seq_dicts[rows[0].strip()] = rows[1].strip()

        for each_ref_acc, list_of_query_acc in grouped_dict.items():
            with open(join(grouped_fasta, each_ref_acc + '.fasta'), 'a') as write_file:
                for each_query_acc in list_of_query_acc:
                    seqs = seq_dicts[each_query_acc]
                    write_file.write(">" + each_query_acc + '\n')
                    for i in range(0, len(seqs), 80):
                        write_file.write(seqs[i:i + 80] + '\n')

        ref_seqs = fasta(db_path)
        for rows in ref_seqs:
            accs, accs_type = rows[0].split("|")
            seq_dicts[accs] = rows[1].strip()
            if accs_type == "master":
                master_seq_dict[accs] = rows[1].strip()

        for ref_acc_info in grouped_dict.keys():
            if "|" in ref_acc_info:
                split_info = ref_acc_info.split("|")
                ref_accs = split_info[0]
                acc_type = split_info[1]
                seqs = seq_dicts[ref_accs]
            else:
                ref_accs = ref_acc_info
                acc_type = ""
                seqs = seq_dicts[ref_acc_info]

            with open(join(ref_seq_dir, ref_accs + '.fasta'), 'w') as write_file:
                write_file.write(">" + ref_accs + '\n')
                for i in range(0, len(seqs), 80):
                    write_file.write(seqs[i:i + 80] + '\n')

        #writing master sequences
        master_fasta = ""
        for master_acc, seqs in master_seq_dict.items():
            master_fasta = join(master_seq_dir, master_acc + '.fasta')
            with open(join(master_seq_dir, master_acc + '.fasta'), 'w') as write_file:
                write_file.write(">" + master_acc + '\n')
                for i in range(0, len(seqs), 80):
                    write_file.write(seqs[i:i + 80] + '\n')


        #merge reference sequence
        write_file = open(ref_path, 'w')
        for each_ref_fa in os.listdir(ref_seq_dir):
            df = fasta(join(ref_seq_dir, each_ref_fa))
            for rows in df:
                write_file.write(">" + rows[0].strip() + "\n" + rows[1].strip() + "\n")
        write_file.close()


        for each_ref in os.listdir(ref_seq_dir):
            os.system('cat ' + master_fasta + ' ' + f'{join(ref_seq_dir, each_ref)}' + ' >' f'{join(master_and_reference_merged, each_ref)}') 
            print(f"Merging complete {join(master_and_reference_merged, each_ref)}") 


    def blastn(self, tmp_dir, query_path, db_path):
        db_file_name = os.path.basename(db_path)
        output_file = join(tmp_dir, "query_tophits.tsv")
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

    def mafft_ref_sequences(self, ref_acc_file, output_dir):
        output_file = self.path_to_basename(ref_acc_file) + "_aln.fasta"
        output_path = join(output_dir, output_file)
        try:
            with open(output_path, 'w') as out_f:
                subprocess.run(['mafft', ref_acc_file], stdout=out_f, check=True)
                print(f"mafft ran successfully on {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error running MAFFT: {e}")

    def mafft_query_sequences(self, query_file, ref_alignment_file, output_dir):
        output_file = self.path_to_basename(ref_alignment_file) + "_aln.fasta"
        output_path = join(output_dir, output_file)
        ref_alignment_file = ref_alignment_file.split('.')[0] + "_aln.fasta"
        try:
            with open(output_path, 'w') as out_f:
                subprocess.run(['mafft', '--add', query_file, '--keeplength', ref_alignment_file], stdout=out_f, check=True)
                print(f"mafft ran successfully on {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error running MAFFT: {e}")

    def extract_matching_sequences(self, fasta_dir, output_file=None):
        results = {}

        for file in os.listdir(fasta_dir):
            if file.endswith("_aln.fasta"):
                accession = file.split("_")[0]
                file_path = os.path.join(fasta_dir, file)

                for record in SeqIO.parse(file_path, "fasta"):
                    if accession in record.id:
                        results[accession] = str(record.seq)
                        break

        if output_file:
            with open(output_file, "w") as out_f:
                for acc, seq in results.items():
                    out_f.write(f">{acc}\n{seq}\n")

    def __extract_ref_seq(self, db_path):

        db_path = join("analysis")
        os.makedirs(db_path, exist_ok=True)
        write_file = open(join(db_path, "db.fa"), 'w')

        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT primary_accession, accession_type FROM meta_data where accession_type='reference' or accession_type='master'")
            ref_accs = cursor.fetchall()
            ref_accs_list = [(item[0], item[1]) for item in ref_accs]

            for each_acc, accession_type in ref_accs_list:
                cursor.execute("SELECT sequence FROM sequences WHERE header = %s", [each_acc])
                if each_acc == "NC_001542": accession_type="master"
                result = cursor.fetchone()
                if result:
                    sequence = result[0]
                    acc_type = "|" + accession_type
                    write_file.write(">" + each_acc + acc_type)
                    write_file.write("\n")
                    write_file.write(sequence)
                    write_file.write("\n")
            write_file.close()

        return 
    
    def __run_makeblastdb(self, db_fasta):
        db_file_name = db_fasta
        command = [
            'makeblastdb',
            '-in', db_fasta,
            '-out', join(db_file_name),
            '-title', "alignment",
            '-dbtype', 'nucl'
        ]
        try:
            subprocess.run(command, check=True)
            print(f"makeblastdb ran successfully on {db_file_name}")
        except subprocess.CalledProcessError as e:
            print(f"Error running makeblastdb: {e}")


if __name__ == "__main__":

	processor = Tasks(database="RABV_NEW")
	processor.process_genbank_submission()
     


    # def process_genbank_submission(self):
    #     db_path = Path(join("analysis","db.fa"))
    #     print(db_path)

    #     job = "test"
    #     os.makedirs(join("jobs",job), exist_ok=True)

    #     # write_fasta_file(sequences, "jobs/"+job+"/input.fasta")


    #     tmp_dir = "jobs/"+job

    #     if not db_path.is_file():
    #         self.__extract_ref_seq()

    #     # self.__blast_analysis(tmp_dir)
    #     # self.__next_align_analysis(tmp_dir)
    #     # self.__pad_alignment(tmp_dir)
    #     # self.__query_feature_table(tmp_dir)

    #     return



    # # PRIVATE FUNCTIONS


    # def __blast_analysis(self, tmp_dir):


    #     values = {}
    #     records = {}

    #     query_path = join(tmp_dir, "input.fasta")
    #     db_path = join("analysis", "db.fa")
    #     query_tophits = join(tmp_dir, "query_tophits.tsv")
    #     query_tophits_uniq = join(tmp_dir, "query_tophits_uniq.tsv")
    #     sorted_fasta = join(tmp_dir, "sorted_fasta")
    #     merged_fasta = join(tmp_dir, "merged_fasta")
    #     sorted_all = join(tmp_dir, "sorted_all")
    #     grouped_fasta = join(tmp_dir, "grouped_fasta")
    #     ref_seq_dir = join(tmp_dir, "reference_sequences")
    #     os.makedirs(sorted_fasta)
    #     os.makedirs(merged_fasta)
    #     os.makedirs(sorted_all)
    #     os.makedirs(grouped_fasta)
    #     os.makedirs(ref_seq_dir)

    #     self.__run_makeblastdb(db_path)
    #     self.__blastn(tmp_dir, query_path, db_path)

    #     with open(query_tophits, newline='') as file:
    #         reader = csv.reader(file, delimiter='\t')
    #         for row in reader:
    #             col1, col2, col3, col4 = row[0], row[1], float(row[2]), row[3]
    #             if col1 in values:
    #                 existing_value = values[col1]
    #                 if col3 > existing_value:
    #                     records[col1] = [col1, col2, col3, col4]
    #                     values[col1] = col3
    #             else:
    #                 records[col1] = [col1, col2, col3, col4]
    #                 values[col1] = col3

    #     with open(query_tophits_uniq, 'w', newline='') as file:
    #         writer = csv.writer(file, delimiter='\t')
    #         for record in records.values():
    #             writer.writerow(record)

    #     seq_dicts = {}
    #     query_seqs = fasta(query_path)

    #     for rows in query_seqs:
    #         seq_dicts[rows[0].strip()] = rows[1].strip()

    #     #Seperate plus and minus strand sequences
    #     for each_line in open(query_tophits_uniq, 'r'):
    #         query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
    #         if strand == "plus":
    #             with open(join(sorted_fasta, "plus.fa"), 'a') as file_plus:
    #                 file_plus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")
    #         else:
    #             with open(join(sorted_fasta, "minus.fa"), 'a') as file_minus:
    #                 file_minus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")

    #     for each_file in os.listdir(sorted_fasta):
    #         if "minus" in each_file:
    #             command = ["seqkit", "seq", "-r", "-p", "-v", "-t", "dna", join(sorted_fasta, each_file), ">", join(merged_fasta, each_file)]
    #         else:
    #             command = ["cp", join(sorted_fasta, each_file), join(merged_fasta, each_file)]

    #         try:
    #             print(' '.join(command))
    #             os.system(' '.join(command))
    #             print(f"seqkit ran successfully for {each_file}")
    #         except subprocess.CalledProcessError as e:
    #             print(f"Error running seqkit: {e}")

    #     file_list = []
    #     for each_file in os.listdir(merged_fasta):
    #         prefix = merged_fasta + "/"
    #         file_list.append(prefix + each_file)

    #     command = ["cat", " ".join(file_list), ">", join(sorted_all, "query_seq.fa")]
    #     print(' '.join(command))
    #     try:
    #         os.system(' '.join(command))
    #         print(f"concatenation sucessful")
    #     except subprocess.CalledProcessError as e:
    #         print(f"Error in concatenation: {e}")

    #     grouped_dict = {}
    #     for each_line in open(query_tophits_uniq, 'r'):
    #         query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
    #         if ref_acc not in grouped_dict:
    #             grouped_dict[ref_acc] = [query_acc]
    #         else:
    #             grouped_dict[ref_acc].append(query_acc)

    #     seq_dicts = {}
    #     query_seqs = fasta(join(sorted_all, "query_seq.fa"))
    #     for rows in query_seqs:
    #         seq_dicts[rows[0].strip()] = rows[1].strip()

    #     for each_ref_acc, list_of_query_acc in grouped_dict.items():
    #         with open(join(grouped_fasta, each_ref_acc + '.fasta'), 'a') as write_file:
    #             for each_query_acc in list_of_query_acc:
    #                 seqs = seq_dicts[each_query_acc]
    #                 write_file.write(">" + each_query_acc + '\n')
    #                 for i in range(0, len(seqs), 80):
    #                     write_file.write(seqs[i:i + 80] + '\n')

    #     ref_seqs = fasta(db_path)
    #     for rows in ref_seqs:
    #         seq_dicts[rows[0].strip()] = rows[1].strip()

    #     for ref_accs in grouped_dict.keys():
    #         seqs = seq_dicts[ref_accs]
    #         with open(join(ref_seq_dir, ref_accs + '.fasta'), 'w') as write_file:
    #             write_file.write(">" + ref_accs + '\n')
    #             for i in range(0, len(seqs), 80):
    #                 write_file.write(seqs[i:i + 80] + '\n')

    #     return 
    	
    # def __next_align_analysis(self, tmp_dir):

    #     query_dir = join(tmp_dir, "grouped_fasta")
    #     ref_dir = join(tmp_dir, "reference_sequences") 
    #     query_aln_output_dir = join(tmp_dir, "query_aln")

    #     for each_query_file in os.listdir(query_dir):
    #         ref_file = each_query_file
    #         self.nextalign(
    #             join(query_dir, each_query_file),
    #             join(ref_dir, ref_file),
    #             query_aln_output_dir
    #         )

    #     return 
    
    # def __pad_alignment(self, tmp_dir):
    #     os.makedirs(join(tmp_dir, "pad_alignment"), exist_ok=True)
    #     nextalign_op_path = join(tmp_dir,"query_aln")
    #     master_seq_len = 11932
    #     query_acc_lst = []
    #     ref_coords = self.load_ref_aln_table()
    #     write_file = open(join(tmp_dir, "pad_alignment", "paded_query.fa"), 'w')
    #     for each_query_alignment in os.listdir(nextalign_op_path):
    #         print(each_query_alignment)
    #         ref_acc = each_query_alignment
    #         fasta_seqs_obj = fasta(join(nextalign_op_path, each_query_alignment, ref_acc + ".aligned.fasta"))
    #         for row in fasta_seqs_obj:
    #             header = row[0]
    #             seq = row[1]
    #             start, end = ref_coords[ref_acc][0][0], ref_coords[ref_acc][0][1]
    #             seq_len = len(seq)
    #             if header not in query_acc_lst:
    #                 query_acc_lst.append(header)
    #                 if len(ref_coords[ref_acc]) <= 1:
    #                     ref_cord_diff = abs(int(ref_coords[ref_acc][0][0]) - int(ref_coords[ref_acc][0][1]))
    #                     if seq_len != ref_cord_diff:
    #                         start, end = ref_coords[ref_acc][0][0], seq_len
    #                     prefix_char = "-" * (int(start)-1)
    #                     suffix_char = "-" * abs(len(prefix_char) + seq_len - master_seq_len)
    #                     write_file.write('>' + header + "\n")
    #                     write_file.write(prefix_char + seq.strip() + suffix_char + "\n")
    #                 else:
    #                     ref_cord_diff = abs(int(ref_coords[ref_acc][0][0]) - int(ref_coords[ref_acc][-1][1]))
    #                     if seq_len != ref_cord_diff + 1:
    #                         start, end = ref_coords[ref_acc][0][0], seq_len
    #                     prefix_char = "-" * (int(start)-1)
    #                     suffix_char = "-" * abs(len(prefix_char) + seq_len - master_seq_len)
    #                     write_file.write('>' + header + "\n")
    #                     write_file.write(prefix_char + seq.strip() + suffix_char + "\n")
    #     write_file.close()

    #     return 
    
    # def __query_feature_table(self, tmp_dir):

    #     sequence_object = fasta(join(tmp_dir, "pad_alignment", "paded_query.fa"))

    #     for row in sequence_object:
    #         content = self.calculate_alignment_coords(row[0], row[1], 30)
    #         for each_cords in content['aligned']:
    #             #write_file = open(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp", row[0] + ".tablel"), 'w')
    #             write_tbl = open(join(tmp_dir, row[0] + ".tbl"), 'w')
    #             start, end = each_cords[0], each_cords[1]
    #             # gff_info = self.find_cds_for_coordinates(gff_dict, start, end)
    #             write_tbl.write(">Feature " + row[0] + "\n")
    #             # for match in gff_info:
    #             # 	cds_start = match['start']
    #             # 	cds_end = match['end']
                    
    #             # 	recalc_cords = self.recalculate_cds_cords(row[1][match['start']-1:match['end']], match['start'], match['end'])
    #             # 	table2asn_cords = self.table2asn_coordinates(gff_dict, recalc_cords[0], recalc_cords[1])
    #             # 	gene = self.extract_gene_name(table2asn_cords[2])
    #             # 	data = [row[0], str(start), str(end), str(recalc_cords[0]), str(recalc_cords[1]), match['product']]					
    #             # 	data_tbl = [str(table2asn_cords[0]), str(table2asn_cords[1])]

    #             # 	write_tbl.write("\t".join(data_tbl) + "\t"  + "gene" + "\n")
    #             # 	write_tbl.write("\t" + "\t" + "\t" + "gene" + "\t" + gene + "\n")
    #             # 	write_tbl.write("\t".join(data_tbl) + "\t"  + "CDS" + "\n")
    #             # 	write_tbl.write("\t" + "\t" + "\t" + "product" + "\t" + table2asn_cords[2] + "\n")
    #             # 	write_tbl.write("\t" + "\t" + "\t" + "gene" + "\t" + gene + "\n")
    #                 #write_file.write('\t'.join(data))
    #                 #write_file.write("\n")
    #     #write_file.close()
    #     write_tbl.close()

    #     return 
    
    # def find_cds_for_coordinates(self, gff_dict, query_start, query_end):
    #     matching_cds = []
    #     for feature, cds_list in gff_dict.items():
    #         if feature=='CDS':
    #             for cds in cds_list:
    #                 cds_start = int(cds['start'])
    #                 cds_end = int(cds['end'])
    #                 if query_end >= cds_start and query_start <= cds_end:
    #                     overlap_start = max(query_start, cds_start)
    #                     overlap_end = min(query_end, cds_end)
    #                     matching_cds.append({
    #                         'start': overlap_start,
    #                         'end': overlap_end,
    #                         'product': cds['product']
    #                 })

    #     return matching_cds
    

    # def __blastn(self, tmp_dir, query_path, db_path):
    #     db_file_name = os.path.basename(db_path)
    #     output_file = join(tmp_dir, "query_tophits.tsv")
    #     command = [
    #         'blastn',
    #         '-query', query_path,
    #         '-db', join("analysis", db_file_name),
    #         '-task', 'blastn',
    #         '-max_target_seqs', '1',
    #         '-max_hsps', '1',
    #         '-out', output_file,
    #         '-outfmt', "6 qacc sacc pident sstrand"
    #     ]
    #     try:
    #         subprocess.run(command, check=True)
    #         print(f"blastn ran successfully. Results saved in {output_file}")
    #     except subprocess.CalledProcessError as e:
    #         print(f"Error running blastn: {e}")
    
    # def nextalign(self, query_acc_path, ref_acc_path, query_aln_op):
    #     accession = self.path_to_basename(query_acc_path)
    #     command = [
    #         'nextalign', 'run',
    #         '--min-seeds', '44',
    #         '--seed-spacing', '50',
    #         '--min-match-rate', '0.1',
    #         '--input-ref', ref_acc_path,
    #         '--output-all', join(query_aln_op, f'{accession}'),
    #         '--output-basename', f'{accession}',
    #         '--include-reference',
    #         query_acc_path
    #     ]

    #     command_str = " ".join(command)
    #     print(f"Executing command: {command_str}")

    #     return_code = os.system(command_str)
    #     if return_code == 0:
    #         print(f"{accession} completed successfully.")
    #     else:
    #         print(f"{accession} failed with return code {return_code}")   

    # def load_ref_aln_table(self):
    #     ref_cords = {}

    #     db_path = join("analysis")
    #     os.makedirs(db_path, exist_ok=True)
    #     write_file = open(join(db_path, "db.fa"), 'w')

    #     with connections[self.database].cursor() as cursor:
    #         cursor.execute("SELECT * FROM features")
    #         reference_features = cursor.fetchall()

    #     for each_item in reference_features:
    #         acc, ref, aln_start, aln_end, cds_start, cds_end, product = each_item
    #         if acc not in ref_cords:
    #             ref_cords[acc] = [(aln_start, aln_end)]
    #         else:
    #             ref_cords[acc].append((aln_start, aln_end))
    #     return ref_cords
    
    # def gap_index(self, sequence):
    #     gap_indices = []
    #     current_gap = []

    #     for index, char in enumerate(sequence, start=1):
    #         if char == '-':
    #             if not current_gap:
    #                 current_gap.append(index)
    #         else:
    #             if current_gap:
    #                 current_gap.append(index - 1)
    #                 gap_indices.append(current_gap)
    #                 current_gap = []

    #     if current_gap:
    #         current_gap.append(len(sequence))
    #         gap_indices.append(current_gap)

    #     return gap_indices

    # def calculate_alignment_coords(self, header, seq, threshold):
    #     gaps = self.gap_index(seq)
    #     seq_len = len(seq)
    #     if len(gaps) == 0:
    #         return {'aligned': [[1, seq_len]], 'unaligned': []}

    #     overall_start = 1
    #     overall_end = gaps[-1][1]

    #     large_gaps = [gap for gap in gaps if (gap[1] - gap[0] + 1) >= threshold]
    #     aligned_segments = []
    #     current_start = overall_start

    #     if len(large_gaps) != 0:
    #         for gap in large_gaps:
    #             gap_start, gap_end = gap
    #             if current_start < gap_start:
    #                 aligned_segments.append([current_start, gap_start - 1])
    #             current_start = gap_end + 1
    #     else:
    #         return {'aligned': [[1, seq_len]], 'unaligned': []}

    #     return {"aligned": aligned_segments, "unaligned": large_gaps}
    