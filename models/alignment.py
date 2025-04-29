from django.db import connections
from models.helper_functions import dictfetchall  # Ensure this function is imported
from models.codon_labeling import get_kuiken2006_codon_labeling
import time
import logging
import django_rq
import rq
import csv

import os
from os.path import join


logger = logging.getLogger(__name__)

# from models.alignment.blast_alignment import BlastAlignment
# from models.NextalignAlignment import NextalignAlignment
# from models.PadAlignment import PadAlignmentSequences
# from models.alignment.GenerateTables import SequenceProcessor
from models.sequences import Sequences


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

class Alignment:
    def __init__(self, database, sequences=None, region=None, nucleotide_or_codon=None, start_coordinate=None, end_coordinate=None, reference_sequence=None):
        self.database = database  

        self.sequences = sequences
        self.region = region
        self.nucleotide_or_codon = nucleotide_or_codon
        self.start_coordinate = start_coordinate
        self.end_coordinate = end_coordinate
        self.reference_sequence = reference_sequence

    def get_reference_sequences_meta_data(self):

        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT * FROM meta_data WHERE primary_accession IN (SELECT alignment_name FROM sequence_alignment)")

            references = dictfetchall(cursor)

        return references

    def get_reference_sequence(self, primary_accession):

        result = {}

        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT * FROM features WHERE accession = %s", [primary_accession])

            features = dictfetchall(cursor)
        
        if not features:
            raise ValueError("Reference sequence with primary_accession {primary_accession} not found")
        
        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT length FROM meta_data WHERE primary_accession = %s", [primary_accession])
            max_whole_genome = cursor.fetchone()[0]

        
        for feature in features:
            codons = get_kuiken2006_codon_labeling(feature["cds_start"], feature["cds_end"])
            feature["codon_start"] = codons[0]
            feature["codon_end"] = codons[1]
            
        
        # features.append({'ref_seq_name':primary_accession, 'feature_name':None, "ref_start":1, "ref_end":max_whole_genome, "product":"Whole Genome", "protein_id":None})
        result["features"] = features

        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT * FROM genes")
            result["genes"] = dictfetchall(cursor)

            cursor.execute("SELECT * FROM sequence_alignment WHERE alignment_name = %s", [primary_accession])
            result["aligned_sequences"] = dictfetchall(cursor)

            cursor.execute("SELECT sequence FROM sequences WHERE header=%s", [primary_accession])
            result["sequence"] = dictfetchall(cursor)[0]

        return result


    #Needs to add error is stuff is missing
    def get_alignments(self):

        if not self.reference_sequence:
            self.reference_sequence = 'NC_001542'
        
        alignments = self.__get_alignments()

        master_alignment = self.__get_master_alignment()

        if not master_alignment:
            raise ValueError("Reference alignment can not be found")
        
        results = self.__parse_alignments(alignments, master_alignment)

        return results

    def build_alignment_fasta_file(self, alignments):

        ofile = open("my_fasta.fasta", "w")

        for alignment in alignments:
            ofile.write(alignment)
        ofile.close()

        return
    



    # def runGenBankSubmitter(self):

    #     parser = ArgumentParser(description='Performs the nextalign of each sequence')
	# parser.add_argument('-q', '--sequence_dir', help='Sequence file directory, it can be single or multiple fasta sequencce files.', required=True)
	# parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp")
	# parser.add_argument('-o', '--output_dir', help='Output directory where processed data and results are stored', default='Table2asn')
	# parser.add_argument('-m', '--metadata', help='A tab limited meta data file name', required=True)
	# parser.add_argument('-n', '--ncbi_submission_template', help='NCBI submission template file generated from "https://submit.ncbi.nlm.nih.gov/genbank/template/submission/"', required=True)
	# parser.add_argument('-gff', '--gff_file', help='Master reference GFF3 file', required=True)
	# parser.add_argument('-db', '--vgtk-db', help='VGTK Database', required=True)

	
	# sequence_dir = 
    #     tmp_dir = 
    #     output_dir = 
    #     metadata = 
    #     ncbi_submission_template = 
    #     gff_file = 
    #     vgtk_db = 
    #     gaps_to_ignore = 30

	# processor = GenBankSequenceSubmitter(args.sequence_dir, args.tmp_dir, args.output_dir, args.metadata, args.ncbi_submission_template, args.gff_file, args.vgtk_db, args.gaps_to_ignore)
	# processor.process()


    #     os.makedirs(join(self.tmp_dir, 'analysis'), exist_ok=True)
	# 	os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'sorted_fasta'), exist_ok=True)
	# 	os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'merged_fasta'), exist_ok=True)
	# 	os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'sorted_all'), exist_ok=True)
	# 	os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'grouped_fasta'), exist_ok=True)
	# 	os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'reference_sequences'), exist_ok=True)
	# 	query_aln_output_dir = join(self.tmp_dir, self.analysis_dir, "query_aln")
	# 	os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'reference_sequence'), exist_ok=True)
	# 	os.makedirs(join(self.tmp_dir, "analysis", self.analysis_dir, "pad_alignment"), exist_ok=True)
		
	# 	self.data_prep()
	# 	self.extract_ref_seq()
	# 	self.blast_analysis()
	# 	self.nextalign_analysis()	
	# 	self.pad_alignment()
	# 	self.query_feature_table()
	# 	self.run_table2asn()
	# 	self.move_sqn_to_output_dir()
	# 	'''
	# 	1. Read all the sequences from the directory (Done)
	# 	2. Read the meta data (Done)
	# 		1. Check if the sequence header is matching to the meta data provided (Done)
	# 	3. Prepare the meta data file suitable to table2asn (Done but needs minor change in header) 
	# 	4. Add them to tmp dir (Done)
	# 		1. Create column seq_type in gB_matrix (Done)
	# 		2. Seperate sequences from gB_matrix (Done)
	# 		3. Method to extract only reference sequeces from DB using sqlite (Done)
	# 	5. Blast the in-house sequences against the reference sequence (Done)
	# 	6. Nextalign the reference and in-house sequences (Done)
	# 	7. Generate feature coordinates for query sequences(Done)
	# 	8. perform table2asn on the tmp directory (Done)
	# 	'''




    def process_tsv_and_add_metadata(self, tsv_path):
        updated_rows = []

        # Open the existing TSV file
        with open(tsv_path+"/query_uniq_tophits.tsv", "r") as file:
            reader = csv.reader(file, delimiter="\t")
            for row in reader:
                sequence_id = row[1]  # Extract the second column (DB search key)
                sequences_helper = Sequences(database='RABV_NEW', primary_accession=sequence_id)
                data = sequences_helper.get_sequence_meta_data()
                # Append metadata to the row
                row.append(data["meta_data"]["length"])
                updated_rows.append(row)

        with open(tsv_path+"/blast_results.tsv", "w", newline="") as file:
            writer = csv.writer(file, delimiter="\t")
            writer.writerows(updated_rows)

        return 




    # redis-server IN THE COMMAND LINE
    # python manage.py rqworker default MUST RUN THIS TO GET IT TO WORK
    def runSequenceAlignment(self):
        print("THIS SHOULD NOT BE STARTING")
        job = rq.get_current_job()
        print(job.id)


        os.makedirs(join("jobs",job.id), exist_ok=True)
        write_fasta_file(self.sequences, "jobs/"+job.id+"/input.fasta")
        # ofile = open("jobs/"+job.id+"/input.fasta", "w")
        # ofile.write(self.sequences)
        # ofile.close()

        """A simple task that sleeps for n seconds."""
        logger.info(f"Starting blast search")
        


        job.meta['status'] = "running blast"
        job.save_meta()
        query = "jobs/"+job.id+"/input.fasta"
        ref_seq = "models/ef437215.fa"
        
        tmp_folder = "jobs/"+job.id
        output_file = "query_tophits.tsv"
        master_acc = "EF437215"

        try:
            job.meta['status'] = "making blast db"
            job.save_meta()
            processor = BlastAlignment(query, ref_seq, tmp_folder, output_file, master_acc)
            processor.run_makeblastdb()
            processor.run_blastn()
            job.meta['status'] = "running blastn"
            job.save_meta()
            processor.process_non_segmented_virus(tmp_folder+"/"+output_file)
            job.meta['status'] = "processing virus"
        except:
            print("blast has failed")
            job.meta['status'] = "FAILED"
            job.save_meta()
            return "FAILED"


# -- sturcutral mutations --> look into this
        
        
        # job.save_meta()

        self.process_tsv_and_add_metadata(tmp_folder)




        time.sleep(3)
        job.meta['status'] = "blast completed"
        job.save_meta()
        job.meta['blast_file'] = "/Users/dana/CVR/backend/jobs/"+job.id+"query_uniq_tophits.tsv"
        logger.info(f"finished blast")
        job.save_meta()
        logger.info(job.meta)
        time.sleep(3)
        job.meta['status'] = "running next align"

        logger.info(f"running next align")
        job.save_meta()
        time.sleep(3)

        os.makedirs(join("jobs",job.id, "Nextalign"), exist_ok=True)

        # query_dir = tmp_folder+'/grouped_fasta'
        # ref_dir = tmp_folder+'/ref_seqs'
        # ref_fa_file = "ef437215.fa"
        # master_seq_dir = tmp_folder+'/master_seq'
        # tmp_dir = tmp_folder+'/Nextalign'
        # master_ref = "EF437215"
        # table_dir = tmp_folder+"/Tables"
        # gaps_to_ignore = 30

    #     	parser = ArgumentParser(description='Performs the nextalign of each sequence')
	# parser.add_argument('-q', '--query_dir', help='Query file directory.', default="tmp/Blast/grouped_fasta")
	# parser.add_argument('-r', '--ref_dir', help='Reference fasta directory', default="tmp/Blast/ref_seqs")
	# parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp/Nextalign")
	# args = parser.parse_args()
        
        query_dir = tmp_folder+'/grouped_fasta'
        ref_dir = tmp_folder+'/ref_seqs'
        tmp_dir = tmp_folder+'/Nextalign'

        processor = NextalignAlignment(query_dir, ref_dir, tmp_dir)
        processor.process()

        # processor = NextalignAlignment(query_dir, ref_dir, tmp_dir, ref_fa_file, master_seq_dir, master_ref, table_dir, gaps_to_ignore)
        # processor.process()

        logger.info(f"finished next align")
        job.meta['status'] = "NextAlign completed"
        job.save_meta()

        job.meta['status'] = "Padding Alignment"
        job.save_meta()


        os.makedirs(join("jobs",job.id, "Pad-Alignment"), exist_ok=True)
        reference_sequence = "models/ef437215.fa"
        input_dir = tmp_folder+'/Nextalign'
        tmp_dir =  tmp_folder+'/Pad-Alignment'
        keep_intermediate_files = "N"
        output_file = "paded-alignment.fa"

        pad_aln = PadAlignmentSequences(reference_sequence, input_dir, tmp_dir, keep_intermediate_files, output_file)
        pad_aln.process()



        job.meta['status'] = "Generating Features"
        job.save_meta()

        os.makedirs(join("jobs",job.id, "GenBank-matrix"), exist_ok=True)

        
        os.makedirs(join("jobs",job.id, "Tables"), exist_ok=True)
        genbank_matrix = 'models/gB_matrix_raw.tsv'
        tmp_dir = tmp_folder+ '/Tables/'
        blast_hits = tmp_folder+ '/query_uniq_tophits.tsv'
        pad_aln = tmp_folder+ '/Pad-Alignment/paded-alignment.fa'
        processor = SequenceProcessor(genbank_matrix, tmp_dir, blast_hits, pad_aln)
        processor.process()




















        time.sleep(3)
        job.meta['status'] = "done"
        job.save_meta()
        return "WOOOO WE ARE DONE!"










    # PRIVATE FUNCTIONS
     
    def __parse_alignments(self, alignments, master_alignment):

        results= []
        if not self.start_coordinate and not self.end_coordinate:  # Use the full region coordinates

            ref_start = master_alignment[0]["ref_start"]
            ref_end = master_alignment[0]["ref_end"]

            if self.nucleotide_or_codon == "codon":  # User wants a codon
                codon_start, codon_end = get_kuiken2006_codon_labeling(ref_start, ref_end)
        else:  # User chooses the coordinates themselves

            ref_start = int(self.start_coordinate)
            ref_end = int(self.end_coordinate)

            if self.nucleotide_or_codon == "codon":  # User wants a codon
                ref_start = master_alignment[0]["ref_start"]
                ref_end = master_alignment[0]["ref_end"]
                codon_start = int(self.start_coordinate)
                codon_end = int(self.end_coordinate)

        for i in range(len(alignments)):
            sub_seq = alignments[i]["alignment"][ref_start:ref_end+1]

            if (self.nucleotide_or_codon == "codon"):
                codons = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]
                selected_codons = codons[codon_start-1:codon_end]
                sub_seq = ''.join(selected_codons)

            if set(sub_seq) != {"-"}:
                results.append(f">" + alignments[i]["sequence_id"] + "\n" + sub_seq + "\n")
        return results


    def __get_alignments(self):

        placeholders = ', '.join(['%s'] * len(self.sequences))  
        query = f'SELECT s.* FROM sequence s WHERE s.sequence_id IN ({placeholders});'

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, self.sequences)
            alignments = dictfetchall(cursor)

        return alignments
    
    def __get_master_alignment(self):
        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT ref_start, ref_end FROM features WHERE ref_seq_name=%s AND feature_name='CDS' AND product=%s", [self.reference_sequence, self.region])
            master_alignment = dictfetchall(cursor)

        return master_alignment
    
    
