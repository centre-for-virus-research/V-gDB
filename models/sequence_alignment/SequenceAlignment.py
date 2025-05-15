#python GenBankSequenceSubmitter-1.py -q ./../test_data/example-2/fasta-seq/ -m ./../test_data/example-2/metadata.tsv -n ./../test_data/example-2/template.sbt -gff NC_001542.gff3 -db ./../../test_version/TING/tmp-rabv/SqliteDB/rabv_apr0825.db

import os
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
import models.sequence_alignment.read_file as read_file
#from PadAlignment import PadAlignment
from models.sequence_alignment.GffToDictionary import GffDictionary
#from FeatureCalculator  import FeatureCordCalculator 
from django.db import connections
from models.helpers import *
class SequenceAlignment:
	def __init__(self, sequence_dir, tmp_dir, output_dir, db, gaps_to_ignore):
		self.sequence_dir = sequence_dir
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.db = db
		self.gaps_to_ignore = gaps_to_ignore


	@staticmethod
	def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]

	def table2asn(self):
		os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)
		command = [
			'table2asn',
			'-indir', join(self.tmp_dir, self.output_dir, "tmp"),
			'-t', self.ncbi_template
		]
		try:
			subprocess.run(command, check=True)
			print(f"Table2asn ran successfully on {self.squence_dir}")
		except subprocess.CalledProcessError as e:
			print(f"Error running Table2asn: {e}")

	def run_makeblastdb(self, db_fasta):
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

	def blastn(self, query_path, db_path):
		db_file_name = os.path.basename(db_path)
		output_file = join(self.tmp_dir, "analysis", "query_tophits.tsv")
		command = [
			'blastn',
			'-query', query_path,
			'-db', join(self.tmp_dir, "analysis", "DB", db_file_name),
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

	def metadata_header(self):
		headers = ["SeqID", "organism", "genotype", "host", "isolate", "isolation_source", "collection-date", "note", "country", "note"]
		#headers = ["FastaID", "sample_id", "isolation_source", "country" ,"place_sampled" ,"host_species" ,"collection_year","notes", "collection_date"]
		return headers

	def load_metadata(self):
		headers = self.metadata_header() 
		df = {}
		with open(self.metadata, mode="r", encoding="utf-8") as file:
			reader = csv.DictReader(file, delimiter="\t")

			missing_headers = [col for col in headers if col not in reader.fieldnames]
			if missing_headers:
				raise ValueError(f"Missing columns in header: {', '.join(missing_headers)}")

			for row in reader:
				seq_id = row["SeqID"]
				if seq_id in df:
					raise ValueError(f"Duplicate SeqID found: {SeqID}")

				df[seq_id] = {key: row[key] for key in headers if key in row}

		return df

	def load_fasta(self):
		fasta_dict = {}

		for filename in os.listdir(self.sequence_dir):
			if filename.endswith(".fasta") or filename.endswith(".fa"):
				filepath = os.path.join(self.sequence_dir, filename)
				with open(filepath, "r") as file:
					for record in SeqIO.parse(file, "fasta"):
						if record.id in fasta_dict:
							raise ValueError(f"Duplicate accession found: {accession}. Ensure unique accessions.")
						fasta_dict[record.id] = str(record.seq)

		return fasta_dict	

	def extract_ref_seq(self):
		db_path = join(self.tmp_dir, "analysis", "DB")
		print(db_path)
		os.makedirs(db_path, exist_ok=True)
		write_file = open(join(db_path, "db.fa"), 'w')
		# conn = sqlite3.connect(self.db)
		# cursor = conn.cursor()
		with connections[self.db].cursor() as cursor:
			cursor.execute("SELECT primary_accession, accession_type FROM meta_data where accession_type='reference' OR accession_type = 'master'")
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
			
            
				
            
                

	def blast_analysis(self):

		values = {}
		records = {}
        
		query_path = join(self.tmp_dir, "input.fasta")
		print(self.tmp_dir)
		ref_path = join(self.tmp_dir, "analysis", "ref.fa")
		db_path = join(self.tmp_dir, "analysis", "DB", "db.fa")
		query_tophits = join(self.tmp_dir, "analysis", "query_tophits.tsv")
		query_tophits_uniq = join(self.tmp_dir, "analysis", "query_tophits_uniq.tsv")
		sorted_fasta = join(self.tmp_dir, "analysis", "sorted_fasta")
		merged_fasta = join(self.tmp_dir, "analysis", "merged_fasta")
		sorted_all = join(self.tmp_dir, "analysis", "sorted_all")
		grouped_fasta = join(self.tmp_dir, "analysis",  "grouped_fasta")
		ref_seq_dir = join(self.tmp_dir, "analysis", "reference_sequences")
		master_seq_dir = join(self.tmp_dir, 'analysis', 'master_sequences')
		master_and_reference_merged = join(self.tmp_dir, 'analysis', 'master_and_reference_merged')

		self.run_makeblastdb(db_path)
		self.blastn(query_path, db_path)

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
		query_seqs = read_file.fasta(query_path)
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
		query_seqs = read_file.fasta(join(sorted_all, "query_seq.fa"))
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		for each_ref_acc, list_of_query_acc in grouped_dict.items():
			with open(join(grouped_fasta, each_ref_acc + '.fasta'), 'a') as write_file:
				for each_query_acc in list_of_query_acc:
					seqs = seq_dicts[each_query_acc]
					write_file.write(">" + each_query_acc + '\n')
					for i in range(0, len(seqs), 80):
						write_file.write(seqs[i:i + 80] + '\n')
	
		ref_seqs = read_file.fasta(db_path)
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
			df = read_file.fasta(join(ref_seq_dir, each_ref_fa))
			for rows in df:
				write_file.write(">" + rows[0].strip() + "\n" + rows[1].strip() + "\n")
		write_file.close()

	
		for each_ref in os.listdir(ref_seq_dir):
			os.system('cat ' + master_fasta + ' ' + f'{join(ref_seq_dir, each_ref)}' + ' >' f'{join(master_and_reference_merged, each_ref)}') 
			print(f"Merging complete {join(master_and_reference_merged, each_ref)}") 

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

       
	def count_gaps_before_position(self, gap_ranges, position):
		"""Count how many positions are removed before a given alignment position."""
		count = 0
		for start, end in gap_ranges:
			if end < position:
				count += (end - start + 1)
			elif start <= position <= end:
				count += (position - start + 1)
		return count
	def get_gap_ranges(self, sequence):
		gap_ranges = []
		start = None

		for i, char in enumerate(sequence):
			if char == '-':
				if start is None:
					start = i + 1  # Convert to 1-based indexing
			else:
				if start is not None:
					gap_ranges.append([start, i])
					start = None

		if start is not None:
			gap_ranges.append([start, len(sequence)])

		return gap_ranges
	def recalculate_cds_coordinates(self, sequence_id, gap_ranges, cds_list, start_offset):
		adjusted_coords = []
		for cds in cds_list:
			cds_start = int(cds['start'])
			cds_end = int(cds['end'])

			gaps_before_start = self.count_gaps_before_position(gap_ranges, cds_start)
			gaps_before_end = self.count_gaps_before_position(gap_ranges, cds_end)

			adj_start = cds_start - gaps_before_start
			adj_end = cds_end - gaps_before_end

			# Apply start offset (e.g., 71 if first gap ends at 70)
			#adj_start += start_offset
			#adj_end += start_offset
			adj_start = cds_start - gaps_before_start + (start_offset - 1)
			adj_end = cds_end - gaps_before_end + (start_offset - 1)


			adjusted_coords.append([adj_start, adj_end])
		return adjusted_coords
	def find_gaps_in_fasta(self, fasta_file, gff_file):
		gff_dict = GffDictionary(gff_file).gff_dict
		cds_list = gff_dict['CDS']

		data = []

		unique_acc_list = []
		meta_data = self.load_metadata()

		for record in SeqIO.parse(fasta_file, "fasta"):
			if record.id in meta_data:
				if record.id not in unique_acc_list:
					unique_acc_list.append(record.id)

					sequence = str(record.seq)
					gaps = self.get_gap_ranges(sequence)

					# Calculate start offset: just after the first gap
					if gaps and gaps[0][0] == 1:
						start_offset = gaps[0][1] + 1
					else:
						start_offset = 1

					adjusted = self.recalculate_cds_coordinates(record.id, gaps, cds_list, start_offset)

					#print(f">{record.id}")
					#print(adjusted)
					data.append([record.id, adjusted, sequence])
		return data
    

def process(self):
		os.makedirs(join(self.tmp_dir, 'analysis'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'sorted_fasta'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'merged_fasta'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'sorted_all'), exist_ok=True)
		grouped_fasta = join(self.tmp_dir, 'analysis', self.analysis_dir, 'grouped_fasta')
		os.makedirs(grouped_fasta, exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'reference_sequences'), exist_ok=True)
		query_aln_output_dir = join(self.tmp_dir, self.analysis_dir, "query_aln")
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'master_sequences'), exist_ok=True)
		master_and_reference_merged = join(self.tmp_dir, "analysis", self.analysis_dir, "master_and_reference_merged")
		mafft_reference_alignment = join(self.tmp_dir, "analysis", self.analysis_dir, "mafft_reference_alignment")
		reference_alignments =  join(self.tmp_dir, "analysis", self.analysis_dir, "reference_alignments")
		query_ref_alignment = join(self.tmp_dir, "analysis", self.analysis_dir, "query_ref_alignment")
		table2asn_tmp = join(self.tmp_dir, "analysis", self.analysis_dir, "Table2asn", "tmp")
		create_tmp_dir = os.makedirs(table2asn_tmp, exist_ok=True)

		os.makedirs(master_and_reference_merged, exist_ok=True)
		os.makedirs(mafft_reference_alignment, exist_ok=True)
		os.makedirs(reference_alignments, exist_ok=True)
		os.makedirs(query_ref_alignment, exist_ok=True)
					
		self.extract_ref_seq()
		self.blast_analysis() 

		for each_ref_merged_acc in os.listdir(master_and_reference_merged):
			self.mafft_ref_sequences(join(master_and_reference_merged, each_ref_merged_acc), mafft_reference_alignment)

		for each_ref_aln in os.listdir(mafft_reference_alignment):
			output_file = self.path_to_basename(each_ref_aln)
			output_file = join(reference_alignments, output_file + '.fasta')
			self.extract_matching_sequences(join(mafft_reference_alignment), output_file)

		for each_query_seq in os.listdir(grouped_fasta):
			reference_file = self.path_to_basename(each_query_seq)
			self.mafft_query_sequences(join(grouped_fasta, each_query_seq), join(reference_alignments, reference_file), query_ref_alignment)
	
if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the nextalign of each sequence')
	parser.add_argument('-q', '--sequence_dir', help='Sequence file directory, it can be single or multiple fasta sequencce files.', required=True)
	parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp")
	parser.add_argument('-o', '--output_dir', help='Output directory where processed data and results are stored', default='Table2asn')
	parser.add_argument('-m', '--metadata', help='A tab limited meta data file name', required=True)
	parser.add_argument('-n', '--ncbi_submission_template', help='NCBI submission template file generated from "https://submit.ncbi.nlm.nih.gov/genbank/template/submission/"', required=True)
	parser.add_argument('-gff', '--gff_file', help='Master reference GFF3 file', required=True)
	parser.add_argument('-db', '--vgtk-db', help='VGTK Database', required=True)
	parser.add_argument('-gp', '--gaps_to_ignore', help='Lenght of gaps to ignore', default=30)
	
	args = parser.parse_args()

	processor = GenBankSequenceSubmitter(args.sequence_dir, args.tmp_dir, args.output_dir, args.metadata, args.ncbi_submission_template, args.gff_file, args.vgtk_db, args.gaps_to_ignore)
	processor.process()
	