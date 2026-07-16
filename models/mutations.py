from django.db import connections
from models.helpers import * 
from Bio.Seq import Seq
from models.filter_helper import FilterHelper


REFERENCE_SEQUENCES_MAP = {"RABV":"NC_001542"}

class Mutations:
    """
    A class to analyze genetic mutations by comparing sequence alignments
    to a master reference sequence in a given database.
    """

    def __init__(self, database, reference_sequence=None):
        """
        Initialize the Mutations object.

        Args:
            database (str): Django database alias to use for queries.
            reference_sequence (str, optional): Accession ID of the reference sequence. Defaults to 'NC_001542'.
        """
        self.database = database
        self.reference_sequence = reference_sequence or REFERENCE_SEQUENCES_MAP[self.database]

    # def _add_filters(self, params):

    #     where_clauses = []
    #     where_params = []
    #     columns_picked = []
    #     for key, value in params.items():
    #         columns_picked.append(f"hl.{key} as host")
    #         if isinstance(value, list):
    #             placeholders = ', '.join(['%s'] * len(value))
    #             where_clauses.append(f"hl.{key} IN ({placeholders})")
    #             where_params.extend(value)
    #         else:
    #             where_clauses.append(f"hl.{key} = %s")
    #             where_params.append(value)

    #     where_str = ' AND '.join(where_clauses)
    #     columns_str = ', '.join(columns_picked)
    #     return where_str, where_params, columns_str

    def _parse_and_translate_results(self, row):

        primary_accession = row["primary_accession"]
        host = row.get("host", "")
        alignment = row['alignment']
        product = row['product']
        cds_start = row['cds_start']
        cds_end = row['cds_end']

        if alignment is None:
            return None
        
        if cds_start is None or cds_end is None:
            return None
        sub_seq = alignment[cds_start: cds_end+1]

        if len(sub_seq) < 3:
            return None

        clean_seq = sub_seq.replace("-", "N")

        protein = str(Seq(clean_seq).translate(to_stop=False))

        results = {
                    "primary_accession": primary_accession,
                    "host": host,
                    "protein": protein,
                    "region":product
                }
        
        return results
    
    def _get_query(self, columns_str, where_sql, params):
        query = f"""
                    SELECT md.primary_accession, sa.alignment, f.product, f.cds_start, f.cds_end, {columns_str} 
                    FROM meta_data md
                    JOIN host_lineage hl ON hl.taxa_id = md.host_taxa_id 
                    JOIN sequence_alignment sa ON sa.sequence_id = md.primary_accession 
                    JOIN features f on f.accession = sa.sequence_id
                    {where_sql}
                """
        with connections[self.database].cursor() as cursor:
            cursor.execute(query, params)
            results = dictfetchall(cursor)

        return results
    
    def _get_reference_query(self, master_accession):
        query = f"""
                SELECT sa.sequence_id as primary_accession, sa.alignment, f.product, f.cds_start, f.cds_end 
                FROM sequence_alignment sa 
                JOIN features f ON f.accession = sa.sequence_id
                WHERE sa.sequence_id = %s
            """
        with connections[self.database].cursor() as cursor:
            cursor.execute(query, [master_accession])
            results = dictfetchall(cursor)

        return results

        
    def get_host_adaptations(self, params, master_accession):
        where_clauses = []
        if params:
            filterHelper = FilterHelper(filters=params, database=self.database)
            where_str, filter_params, columns_str = filterHelper.add_filters_host_mutations()
            # where_str, filter_params, columns_str = self._add_filters(params)
            where_clauses.append(where_str)

        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        translated_query_results, translated_ref_results = [], []

        query_rows = self._get_query(columns_str, where_sql, filter_params)
        
        for row in query_rows:

            translated_row = self._parse_and_translate_results(row)
            if translated_row:
                translated_query_results.append(translated_row)

        ref_rows = self._get_reference_query(master_accession)
        
        for row in ref_rows:

            translated_row = self._parse_and_translate_results(row)
            if translated_row:
                translated_ref_results.append(translated_row)

        results = {
                    "reference_protein": translated_ref_results, 
                    "translated_sequences": translated_query_results
                }
        return results




    # def get_adaptive_mutations_chart_RABV(self, params):
    #     where_clauses = []
    #     if params:
    #         where_str, filter_params, columns_str = self._add_filters(params)
    #         where_clauses.append(where_str)

    #     where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""


    #     with connections[self.database].cursor() as cursor:
            

    #         query = f"""
    #                     SELECT md.primary_accession, sa.alignment, f.product, f.cds_start, f.cds_end, {columns_str} 
    #                     FROM meta_data md
    #                     JOIN host_lineage hl ON hl.taxa_id = md.host_taxa_id 
    #                     JOIN sequence_alignment sa ON sa.sequence_id = md.primary_accession 
    #                     JOIN features f on f.accession = sa.sequence_id
    #                     {where_sql}
    #                 """
    #         print(query, filter_params)
    #         cursor.execute(query, filter_params)
    #         rows = dictfetchall(cursor)
            
    #         query_ref = f"""
    #                         SELECT sa.alignment, f.product, f.cds_start, f.cds_end 
    #                         FROM sequence_alignment sa 
    #                         JOIN features f ON f.accession = sa.sequence_id
    #                         WHERE sa.sequence_id = %s
    #                     """
    #         cursor.execute(query_ref, ['NC_001542'])
    #         ref_results = dictfetchall(cursor)
        
    #     translated_ref_results = []
    #     for row in ref_results:

    #         primary_accession = 'NC_001542'
    #         alignment = row['alignment']
    #         product = row['product']
    #         cds_start = row['cds_start']
    #         cds_end = row['cds_end']

    #         if alignment is None:
    #             continue
    #         sub_seq = alignment[cds_start: cds_end+1]

    #         cds_seq = sub_seq.replace("-", "")
    #         if len(cds_seq) < 3:
    #             continue
    #         protein = str(Seq(cds_seq).translate(to_stop=False))
            
    #         translated_ref_results.append({
    #             "primary_accession": primary_accession,
    #             "protein": protein,
    #             "region":product
    #         })


    #     translated_results = []
    #     # print(rows)
    #     for row in rows:
    #         # print(row)
    #         primary_accession = row['primary_accession']
    #         host = row['host']
    #         alignment = row['alignment']
    #         product = row['product']
    #         cds_start = row['cds_start']
    #         cds_end = row['cds_end']

    #         if alignment is None:
    #             continue
    #         print(cds_start, cds_end)
    #         if cds_start is None or cds_end is None:
    #             sub_seq = ''
    #         else:
    #             sub_seq = alignment[cds_start: cds_end+1]
    #         # cds_seq = sub_seq.replace("-", "")
    #         # if len(sub_seq) < 3:
    #         #     continue
    #         # protein = str(Seq(sub_seq).translate(to_stop=False))

    #         codons = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]
    #         protein = ''
    #         for codon in codons:
    #             protein += translateCodon(codon)
                

            
            
    #         translated_results.append({
    #             "primary_accession": primary_accession,
    #             "host": host,
    #             "protein": protein,
    #             "region":product
    #         })
    #     print(translated_results)

    #     results = {"reference_protein": translated_ref_results, "translated_sequences": translated_results}
    #     # print(translated_results)
    #     return results
    

    # THIS ONE IS FOR FLU -- NEEDS ACTUAL WORK FOR IT
    def get_adaptive_mutations_chart(self, segment):
        
        pb2_ref_seq = "ATGGAAAGAATAAAAGAACTAAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGATCATATGGCCATAATCAAGAAATACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCAGACAAGAGAATAATGGAGATGATTCCTGAAAGGAATGAGCAAGGACAAACGCTTTGGAGCAAGACAAATGATGCTGGGTCGGACAGAGTGATGGTGTCTCCCCTAGCTGTAACTTGGTGGAACAGGAATGGGCCGACAACAAGTACAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGTTGAAAGGTTAAAACATGGAACCTTCGGTCCCGTTCATTTCCGAAACCAAGTTAAAATACGTCGCCGGGTGGATATAAACCCGGGCCATGCAGATCTCAGTGCTAAAGAAGCACAAGATGTTATCATGGAGGTCGTTTTCCCAAATGAAGTGGGAGCTAGAATATTGACATCAGAGTCGCAATTGACAATAACAAAAGAGAAGAAAGAAGAGCTCCAGGATTGTAAAATTGCTCCTTTAATGGTGGCATACATGTTGGAAAGAGAACTGGTCCGCAAAACCAGATTTCTACCGGTAGCAGGCGGAACAAGCAGTGTGTACATTGAGGTATTGCATTTGACTCAAGGGACCTGTTGGGAACAGATGTACACTCCCGGCGGAGAAGTAAGAAATGATGATGTTGACCAGAGTTTGATCATCGCTGCCAGAAACATTGTTAGGAGAGCAACAGTATCAGCGGACCCACTGGCATCACTCTTGGAGATGTGTCACAGCACACAAATTGGGGGAATAAGGATGGTGGACATCCTTAGGCAAAACCCAACTGAGGAGCAAGCTGTGGATATATGCAAAGCAGCAATGGGTTTGAGGATCAGTTCATCCTTTAGCTTTGGAGGCTTCACTTTCAAAAGAACAAATGGATCATCCGTCAAGAAGGAAGAGGAAGTGCTTACAGGCAACCTCCAAACATTGAAAATAAAAGTACATGAGGGGTATGAAGAATTCACAATGGTTGGGCGGAGAGCAACAGCTATCCTGAGGAAAGCAACTAGAAGGCTGATTCAGTTGATAGTAAGTGGAAGAGATGAACAATCAATCGCTGAAGCGATCATTGTAGCAATGGTGTTCTCACAGGAGGATTGCATGATAAAGGCAGTCCGAGGCGATCTGAATTTCGTGAACAGAGCAAACCAAAGATTGAACCCCATGCATCAACTCCTGAGGCACTTCCAAAAAGATGCAAAAGTGCTGTTTCAGAACTGGGGAATTGAACCTATTGACAATGTCATGGGGATGATCGGAATATTACCTGACATGACTCCAAGCGCAGAGATGTCACTGAGAGGAGTGAGAGTTAGTAAGATGGGAGTAGATGAATATTCCAGCACGGAGAGAGTGGTGGTGAGTATTGACCGTTTCTTGAGGGTCCGAGATCAGCAGGGGAACGTACTCTTATCTCCTGAAGAGGTTAGTGAAACACAGGGAACAGAGAAGTTGACAATAACATATTCATCCTCAATGATGTGGGAAATCAACGGTCCTGAGTCAGTGCTTGTTAACACTTATCAATGGATCATCAGGAATTGGGAGACTGTAAAGATTCAATGGTCTCAAGATCCCACAATGCTGTACAATAAGATGGAGTTTGAATCGTTCCAATCCTTGGTGCCAAAGGCTGCCAGAAGCCAATATAGTGGATTTGTGAGAACACTATTCCAACAGATGCGTGATGTTTTGGGGACATTTGATACTGTCCAAATAATCAAGCTGCTACCATTTGCAGCAGCCCCACCGGAGCCGAGCAGAATGCAGTTTTCTTCTCTAACTGTGAATGTGAGAGGCTCAGGAATGAGAATACTCGTGAGGGGTAACTCCCCCGTGTTCAACTACAACAAGGCAACCAAAAGGCTTACAGTCCTCGGAAAGGACGCAGGTGCATTAACAGAAGATCCAGACGAGGGAACAGCCGGGGTGGAATCTGCAGTATTGAGGGGATTCCTAATTCTAGGCAGAGAGGACAAAAGATATGGACCCGCATTGAGCATCAATGAACTGAGCAATCTTGCAAAAGGGGAGAAGGCTAATGTATTGATAATGCAAGGAGACGTGGTGTTGGTAATGAAACGGAAACGGGACTTTAGCATACTTACTGACAGCCAGACAGCGACCAAAAGAATTCGGATGGCCATCAAT---TAG"
        print("starting")
        with connections["FLUV"].cursor() as cursor:
            # cursor.execute("""SELECT cm.primary_accession, md.host_taxa_id
            #                 FROM cluster_members cm
            #                 JOIN meta_data md ON md.primary_accession = cm.primary_accession
            #                 WHERE cm.segment = %s""", ['PB2'])
            cursor.execute("""SELECT cm.primary_accession, md.host_scientific_name as host, sa.alignment
                            FROM cluster_members cm
                            JOIN meta_data md ON md.primary_accession = cm.primary_accession
                            JOIN sequence_alignment sa ON sa.sequence_id = cm.primary_accession
                            WHERE cm.segment = %s""", ['PB2'])
            print("done")
            rows = dictfetchall(cursor)

        translated_results = []
        # print(rows)
        for row in rows:
            primary_accession = row['primary_accession']
            host = row['host']
            alignment = row['alignment']
            if alignment is None:
                continue

            # Remove gaps if this is an alignment
            cds_seq = alignment.replace("-", "")
            protein = str(Seq(cds_seq).translate(to_stop=False))
            
            translated_results.append({
                "primary_accession": primary_accession,
                "host": host,
                "protein": protein,
            })


        ref_seq = pb2_ref_seq.replace("-", "")

        
        # print(protein)
        ref_protein = str(Seq(ref_seq).translate(to_stop=False))

        results = {"reference_protein": ref_protein, "translated_sequences": translated_results}
        # print(translated_results)
        return results
