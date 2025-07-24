import logging
import rq
import os
from os.path import join
from pathlib import Path
from models.sequence_alignment.SequenceAlignment import *
from models.sequence_alignment.GenBankSequenceSubmitter import *
from models.sequence_alignment.GffToDictionary import GffDictionary
import urllib.parse

def save_query_to_fasta(encoded_query, filename="output.fasta"):
    # Step 1: URL-decode the query
    

    # Step 2: Save to file
    with open(filename, "w") as fasta_file:
        fasta_file.write(decoded_query)

    print(f"Saved FASTA to {filename}")

logger = logging.getLogger(__name__)

# ALL OF THESE USE redis server to work
class Tasks:
    def __init__(self, database, query, tmp_dir=None):
        self.database = database  
        self.query = query
        self.tmp_dir = tmp_dir

    @staticmethod
    def path_to_basename(file_path):
        path = os.path.basename(file_path)
        return path.split('.')[0]
    
        



    def process_genbank_submission(self, job):
        print("HI")
        



        
        

        

    # OFFICEAL SEQUENCE ALIGNMNET TASK CODE
    def run_sequence_alignment(self):
        # Add in sequence alignment logic here
        logger.info(f"Starting sequence alignment")
        job = rq.get_current_job()
        job.meta['status'] = {}
        job.save_meta()
        job.meta['status']["message"] = 'creating job directories'
        job.save_meta()
        os.makedirs(join("jobs", job.id), exist_ok=True)
        self.tmp_dir = join("jobs", job.id)
        print(self.tmp_dir)

        sequence_dir = join("jobs", job.id, 'input.fasta')
        decoded_query = urllib.parse.unquote(self.query)
        with open(join(sequence_dir), 'a') as file_plus:
            file_plus.write(decoded_query)
                                   

        
        os.makedirs(join(self.tmp_dir, 'analysis'), exist_ok=True)

        os.makedirs(join(self.tmp_dir, 'output'), exist_ok=True)
        output_dir = join(self.tmp_dir, 'output')


        os.makedirs(join(self.tmp_dir, 'analysis', 'sorted_fasta'), exist_ok=True)
        os.makedirs(join(self.tmp_dir, 'analysis', 'merged_fasta'), exist_ok=True)
        os.makedirs(join(self.tmp_dir, 'analysis', 'sorted_all'), exist_ok=True)
        grouped_fasta = join(self.tmp_dir, 'analysis', 'grouped_fasta')
        os.makedirs(grouped_fasta, exist_ok=True)
        os.makedirs(join(self.tmp_dir, 'analysis', 'reference_sequences'), exist_ok=True)
        query_aln_output_dir = join(self.tmp_dir, "query_aln")
        os.makedirs(join(self.tmp_dir, 'analysis', 'master_sequences'), exist_ok=True)
        master_and_reference_merged = join(self.tmp_dir, "analysis", "master_and_reference_merged")
        mafft_reference_alignment = join(self.tmp_dir, "analysis", "mafft_reference_alignment")
        reference_alignments =  join(self.tmp_dir, "analysis", "reference_alignments")
        query_ref_alignment = join(self.tmp_dir, "analysis", "query_ref_alignment")
        table2asn_tmp = join(self.tmp_dir, "analysis", "Table2asn", "tmp")
        create_tmp_dir = os.makedirs(table2asn_tmp, exist_ok=True)

        os.makedirs(master_and_reference_merged, exist_ok=True)
        os.makedirs(mafft_reference_alignment, exist_ok=True)
        os.makedirs(reference_alignments, exist_ok=True)
        os.makedirs(query_ref_alignment, exist_ok=True)


        processor = SequenceAlignment(sequence_dir, 
                            self.tmp_dir, 
                            output_dir, 
                            self.database,
                            30)

        job.meta['status']['message'] = 'extracting reference sequences'
        job.save_meta()
        processor.extract_ref_seq()

        job.meta['status']['blast'] = 'in-progress'
        job.save_meta()
        job.meta['status']['message'] = 'starting blast search'
        job.save_meta()

        processor.blast_analysis()
        
        job.meta['status']['blast'] = 'done'
        job.save_meta()
        job.meta['status']['message'] = 'starting MAFFT alignment'
        job.save_meta()
        job.meta['status']['alignment'] = 'in-progress'
        job.save_meta()

        for each_ref_merged_acc in os.listdir(master_and_reference_merged):
            processor.mafft_ref_sequences(join(master_and_reference_merged, each_ref_merged_acc), mafft_reference_alignment)

        for each_ref_aln in os.listdir(mafft_reference_alignment):
            output_file = self.path_to_basename(each_ref_aln)
            output_file = join(reference_alignments, output_file + '.fasta')
            processor.extract_matching_sequences(join(mafft_reference_alignment), output_file)

        for each_query_seq in os.listdir(grouped_fasta):
            reference_file = self.path_to_basename(each_query_seq)
            processor.mafft_query_sequences(join(grouped_fasta, each_query_seq), join(reference_alignments, reference_file), query_ref_alignment)
        

        gff_file = "/Users/dana/CVR/V-gDB_Projects/backend/V-gDB/models/sequence_alignment/NC_001542.gff3"
        gff_dict = GffDictionary(gff_file).gff_dict

        non_redundant_acc = []
        for each_query_ref_aln in os.listdir(query_ref_alignment):
            partial_gaps = processor.find_gaps_in_fasta(join(query_ref_alignment, each_query_ref_aln), gff_file)
            print(partial_gaps)
            print("HELLO")
            for each_acc in partial_gaps:
                print(each_acc)
            
        
                acc_num = each_acc[0]
                coordinates = each_acc[1]

                with open(join(table2asn_tmp, acc_num + ".tbl"), "w") as out:
                    out.write(f">Feature {acc_num}\n")

                    for each_cords in coordinates:
                        cord_start = each_cords[0]
                        cord_end = each_cords[1]
                        print(each_cords)
                        table2asn_cords = processor.table2asn_coordinates(gff_dict, cord_start, cord_end)
                        print(table2asn_cords)
                        if (table2asn_cords != None):
                            start, end, product = table2asn_cords

                            start_clean = start.lstrip("<>")
                            end_clean = end.lstrip("<>")

                            gene_symbol = product.split()[-1]

                            product_name = product.replace(gene_symbol, "").strip(" ,")

                            out.write(f"{start_clean}\t{end_clean}\tgene\n")
                            out.write(f"\t\t\tgene\t{gene_symbol}\n")

                            out.write(f"{start}\t{end}\tCDS\n")
                            out.write(f"\t\t\tproduct\t{product_name}\n")
                            out.write(f"\t\t\tgene\t{gene_symbol}\n")
        job.meta['status']['alignment'] = 'done'
        job.save_meta()
        job.meta['status']['message'] = 'done'
        job.save_meta()


        # # try: 


        # except:
        #     print("this didn't work")
        #     job.meta['status'] = 'failed'
        #     job.save_meta()


        

        # self.process_genbank_submission(job.id)
        



        # job.meta['status']['blast'] = 'in-progress'


        return job.meta
    
    # OFFICEAL SEQUENCE ALIGNMNET TASK CODE
    def run_sequence_alignment_test(self):

        self.tmp_dir = join("jobs", "test")
        print(self.tmp_dir)

        sequence_dir = join("jobs", "test", 'input.fasta')

                                   

        
        # os.makedirs(join(self.tmp_dir, 'analysis'), exist_ok=True)

        # os.makedirs(join(self.tmp_dir, 'output'), exist_ok=True)
        output_dir = join(self.tmp_dir, 'output')
        analysis_dir = str(uuid.uuid4())

        # os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir,'sorted_fasta'), exist_ok=True)
        # os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir,'merged_fasta'), exist_ok=True)
        # os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir,'sorted_all'), exist_ok=True)
        # grouped_fasta = join(self.tmp_dir, 'analysis', analysis_dir,'grouped_fasta')
        # os.makedirs(grouped_fasta, exist_ok=True)
        # os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir,'reference_sequences'), exist_ok=True)
        # query_aln_output_dir = join(self.tmp_dir, "query_aln")
        # os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir,'master_sequences'), exist_ok=True)
        # master_and_reference_merged = join(self.tmp_dir, "analysis", analysis_dir,"master_and_reference_merged")
        # mafft_reference_alignment = join(self.tmp_dir, "analysis", analysis_dir,"mafft_reference_alignment")
        # reference_alignments =  join(self.tmp_dir, "analysis", analysis_dir,"reference_alignments")
        # query_ref_alignment = join(self.tmp_dir, "analysis", analysis_dir,"query_ref_alignment")
        # table2asn_tmp = join(self.tmp_dir, "analysis", analysis_dir,"Table2asn", "tmp")
        # create_tmp_dir = os.makedirs(table2asn_tmp, exist_ok=True)

        # os.makedirs(master_and_reference_merged, exist_ok=True)
        # os.makedirs(mafft_reference_alignment, exist_ok=True)
        # os.makedirs(reference_alignments, exist_ok=True)
        # os.makedirs(query_ref_alignment, exist_ok=True)
        os.makedirs(join(self.tmp_dir, 'analysis'), exist_ok=True)
        os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir, 'sorted_fasta'), exist_ok=True)
        os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir, 'merged_fasta'), exist_ok=True)
        os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir, 'sorted_all'), exist_ok=True)
        grouped_fasta = join(self.tmp_dir, 'analysis', analysis_dir, 'grouped_fasta')
        os.makedirs(grouped_fasta, exist_ok=True)
        os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir, 'reference_sequences'), exist_ok=True)
        query_aln_output_dir = join(self.tmp_dir, analysis_dir, "query_aln")
        os.makedirs(join(self.tmp_dir, 'analysis', analysis_dir, 'master_sequences'), exist_ok=True)
        master_and_reference_merged = join(self.tmp_dir, "analysis", analysis_dir, "master_and_reference_merged")
        mafft_reference_alignment = join(self.tmp_dir, "analysis", analysis_dir, "mafft_reference_alignment")
        reference_alignments =  join(self.tmp_dir, "analysis", analysis_dir, "reference_alignments")
        query_ref_alignment = join(self.tmp_dir, "analysis", analysis_dir, "query_ref_alignment")
        table2asn_tmp = join(self.tmp_dir, "analysis", analysis_dir, "Table2asn", "tmp")
        create_tmp_dir = os.makedirs(table2asn_tmp, exist_ok=True)
        os.makedirs(master_and_reference_merged, exist_ok=True)
        os.makedirs(mafft_reference_alignment, exist_ok=True)
        os.makedirs(reference_alignments, exist_ok=True)
        os.makedirs(query_ref_alignment, exist_ok=True)
        
        processor = GenBankSequenceSubmitter(sequence_dir, 
                            self.tmp_dir, 
                            output_dir, 
                            None, 
                            None,
                            gff_file = None,
                            db = self.database,
                            gaps_to_ignore = 30,
                            analysis_dir= analysis_dir)



					
        processor.extract_ref_seq()
        processor.blast_analysis() 

        for each_ref_merged_acc in os.listdir(master_and_reference_merged):
            processor.mafft_ref_sequences(join(master_and_reference_merged, each_ref_merged_acc), mafft_reference_alignment)

        for each_ref_aln in os.listdir(mafft_reference_alignment):
            output_file = processor.path_to_basename(each_ref_aln)
            output_file = join(reference_alignments, output_file + '.fasta')
            processor.extract_matching_sequences(join(mafft_reference_alignment), output_file)

        for each_query_seq in os.listdir(grouped_fasta):
            reference_file = processor.path_to_basename(each_query_seq)
            processor.mafft_query_sequences(join(grouped_fasta, each_query_seq), join(reference_alignments, reference_file), query_ref_alignment)


        gff_file = "/Users/dana/CVR/V-gDB_Projects/backend/V-gDB/models/sequence_alignment/NC_001542.gff3"
        gff_dict = GffDictionary(gff_file).gff_dict

        non_redundant_acc = []
        for each_query_ref_aln in os.listdir(query_ref_alignment):
            partial_gaps = processor.find_gaps_in_fasta(join(query_ref_alignment, each_query_ref_aln), gff_file)
            print(partial_gaps)
            print("HELLO")
            for each_acc in partial_gaps:
            
        
                acc_num = each_acc[0]
                coordinates = each_acc[1]

                with open(join(table2asn_tmp, acc_num + ".tbl"), "w") as out:
                    out.write(f">Feature {acc_num}\n")

                    for each_cords in coordinates:
                        cord_start = each_cords[0]
                        cord_end = each_cords[1]
                        table2asn_cords = processor.table2asn_coordinates(gff_dict, cord_start, cord_end)

                        start, end, product = table2asn_cords

                        start_clean = start.lstrip("<>")
                        end_clean = end.lstrip("<>")

                        gene_symbol = product.split()[-1]

                        product_name = product.replace(gene_symbol, "").strip(" ,")

                        out.write(f"{start_clean}\t{end_clean}\tgene\n")
                        out.write(f"\t\t\tgene\t{gene_symbol}\n")

                        out.write(f"{start}\t{end}\tCDS\n")
                        out.write(f"\t\t\tproduct\t{product_name}\n")
                        out.write(f"\t\t\tgene\t{gene_symbol}\n")


        # processor.extract_ref_seq()

        # processor.blast_analysis()


        # for each_ref_merged_acc in os.listdir(master_and_reference_merged):
        #     processor.mafft_ref_sequences(join(master_and_reference_merged, each_ref_merged_acc), mafft_reference_alignment)

        # for each_ref_aln in os.listdir(mafft_reference_alignment):
        #     output_file = self.path_to_basename(each_ref_aln)
        #     output_file = join(reference_alignments, output_file + '.fasta')
        #     processor.extract_matching_sequences(join(mafft_reference_alignment), output_file)

        # for each_query_seq in os.listdir(grouped_fasta):
        #     reference_file = self.path_to_basename(each_query_seq)
        #     processor.mafft_query_sequences(join(grouped_fasta, each_query_seq), join(reference_alignments, reference_file), query_ref_alignment)
        

        return 


