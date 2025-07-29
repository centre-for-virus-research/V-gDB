from django.db import connections
import csv

from models.helpers import *
from models.sequences import Sequences


class Alignment:
    def __init__(self, database, sequences=None, region=None, nucleotide_or_codon=None, start_coordinate=None, end_coordinate=None, reference_sequence=None):
        self.database = database  

        self.sequences = sequences
        self.region = region
        self.nucleotide_or_codon = nucleotide_or_codon
        self.start_coordinate = start_coordinate
        self.end_coordinate = end_coordinate
        self.reference_sequence = reference_sequence

   

    #Needs to add error is stuff is missing
    def get_alignments(self):

        if not self.reference_sequence:
            self.reference_sequence = 'NC_001542'
        
        # alignments = self.__get_alignments()

        # master_alignment = self.__get_master_alignment()

        # if not master_alignment:
        #     raise ValueError("Reference alignment can not be found")
        
        # results = self.__parse_alignments(alignments, master_alignment)

        alignments = self.__get_alignments_and_features()
        results = self.__parse_alignments_new(alignments)

        return results

    def process_tsv_and_add_metadata(self, tsv_path):
        updated_rows = []

        # Open the existing TSV file
        with open(tsv_path+"/query_uniq_tophits.tsv", "r") as file:
            reader = csv.reader(file, delimiter="\t")
            for row in reader:
                sequence_id = row[1]  # Extract the second column (DB search key)
                sequences_helper = Sequences(database='RABV')
                data = sequences_helper.get_sequence_meta_data(sequence_id)
                # Append metadata to the row
                row.append(data["meta_data"]["length"])
                updated_rows.append(row)

        with open(tsv_path+"/blast_results.tsv", "w", newline="") as file:
            writer = csv.writer(file, delimiter="\t")
            writer.writerows(updated_rows)

        return 


    # PRIVATE FUNCTIONS


    def __parse_alignments_new(self, alignments):
        results= []

        for a in alignments:
            if not self.start_coordinate and not self.end_coordinate:  # Use the full region coordinates

                ref_start = a["cds_start"]
                ref_end = a["cds_end"]

                if self.nucleotide_or_codon == "codon":  # User wants a codon
                    codon_start, codon_end = get_codon_labeling(ref_start, ref_end)
            else:  # User chooses the coordinates themselves

                ref_start = int(self.start_coordinate)
                ref_end = int(self.end_coordinate)

                if self.nucleotide_or_codon == "codon":  # User wants a codon
                    ref_start = a["cds_start"]
                    ref_end = a["cds_end"]
                    codon_start = int(self.start_coordinate)
                    codon_end = int(self.end_coordinate)

            sub_seq = a["alignment"][ref_start:ref_end+1]

            if (self.nucleotide_or_codon == "codon"):
                codons = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]
                selected_codons = codons[codon_start-1:codon_end]
                sub_seq = ''.join(selected_codons)

            if set(sub_seq) != {"-"}:
                results.append(f">" + a["sequence_id"] + "\n" + sub_seq + "\n")

        return results



     
    def __parse_alignments(self, alignments, master_alignment):

        results= []
        if not self.start_coordinate and not self.end_coordinate:  # Use the full region coordinates

            ref_start = master_alignment[0]["cds_start"]
            ref_end = master_alignment[0]["cds_end"]

            if self.nucleotide_or_codon == "codon":  # User wants a codon
                codon_start, codon_end = get_codon_labeling(ref_start, ref_end)
        else:  # User chooses the coordinates themselves

            ref_start = int(self.start_coordinate)
            ref_end = int(self.end_coordinate)

            if self.nucleotide_or_codon == "codon":  # User wants a codon
                ref_start = master_alignment[0]["cds_start"]
                ref_end = master_alignment[0]["cds_end"]
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
        query = f'SELECT s.* FROM sequence_alignment s WHERE s.sequence_id IN ({placeholders});'
        with connections[self.database].cursor() as cursor:
            cursor.execute(query, self.sequences)
            alignments = dictfetchall(cursor)

        return alignments
    
    def __get_master_alignment(self):
        with connections[self.database].cursor() as cursor:
            if(self.region != 'entirety'):
                cursor.execute("SELECT cds_start, cds_end FROM features WHERE accession=%s AND product=%s", [self.reference_sequence, self.region])
                master_alignment = dictfetchall(cursor)
            else:
                cursor.execute("SELECT length FROM meta_data WHERE primary_accession=%s", [self.reference_sequence])
                tmp = dictfetchall(cursor)
                master_alignment = [{"cds_start": 1, "cds_end":tmp[0]["length"]}]

        return master_alignment


    def __get_alignments_and_features(self):
        placeholders = ', '.join(['%s'] * len(self.sequences))  
        query = f'SELECT s.*, f.* \
                    FROM sequence_alignment s \
                    JOIN features f on f.accession = s.sequence_id \
                    WHERE s.sequence_id IN ("MT862689", "JX987734", "DQ468335") \
                    AND f.product = "transmembrane glycoprotein G"'
        with connections[self.database].cursor() as cursor:
            cursor.execute(query)
            alignments = dictfetchall(cursor)

        return alignments

        
    