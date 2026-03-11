from django.db import connections
from models.helpers import * 
from Bio.Seq import Seq

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
        self.reference_sequence = reference_sequence or 'NC_001542'

    def _add_filters(self, params):

        where_clauses = []
        where_params = []
        columns_picked = []
        for key, value in params.items():
            columns_picked.append(f"hl.{key} as host")
            if isinstance(value, list):
                placeholders = ', '.join(['%s'] * len(value))
                where_clauses.append(f"hl.{key} IN ({placeholders})")
                where_params.extend(value)
            else:
                where_clauses.append(f"hl.{key} = %s")
                where_params.append(value)
        print(where_clauses, where_params)

        where_str = ' AND '.join(where_clauses)
        columns_str = ', '.join(columns_picked)
        return where_str, where_params, columns_str

    def get_adaptive_mutations(self):

        query = "SELECT * FROM adaptive_mutations;"

        with connections[self.database].cursor() as cursor:
            cursor.execute(query)
            mutations = dictfetchall(cursor)

        return mutations

    def get_adaptive_mutations_chart_RABV(self, params):
        where_clauses = []
        if params:
            where_str, filter_params, columns_str = self._add_filters(params)
            where_clauses.append(where_str)

        # query = 'SELECT DISTINCT(class) FROM host_lineage WHERE  ORDER BY class ASC'
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""
        
        # print(query, filter_params)
        print("COLUMNS", columns_str)
        print("SQL", where_sql)

        with connections[self.database].cursor() as cursor:
            

            query = f"""
                        SELECT md.primary_accession, sa.alignment, f.product, f.cds_start, f.cds_end, {columns_str} 
                        FROM meta_data md
                        JOIN host_lineage hl ON hl.taxa_id = md.host_taxa_id 
                        JOIN sequence_alignment sa ON sa.sequence_id = md.primary_accession 
                        JOIN features f on f.accession = sa.sequence_id
                        {where_sql}
                    """
            print(query, filter_params)
            cursor.execute(query, filter_params)
            rows = dictfetchall(cursor)
            
            query_ref = f"""
                            SELECT sa.alignment, f.product, f.cds_start, f.cds_end 
                            FROM sequence_alignment sa 
                            JOIN features f ON f.accession = sa.sequence_id
                            WHERE sa.sequence_id = %s
                        """
            cursor.execute(query_ref, ['NC_001542'])
            ref_results = dictfetchall(cursor)
        
        translated_ref_results = []
        for row in ref_results:

            primary_accession = 'NC_001542'
            alignment = row['alignment']
            product = row['product']
            cds_start = row['cds_start']
            cds_end = row['cds_end']

            if alignment is None:
                continue
            sub_seq = alignment[int(cds_start): int(cds_end)+1]

            cds_seq = sub_seq.replace("-", "")
            if len(cds_seq) < 3:
                continue
            protein = str(Seq(cds_seq).translate(to_stop=False))
            
            translated_ref_results.append({
                "primary_accession": primary_accession,
                "protein": protein,
                "region":product
            })


        translated_results = []
        # print(rows)
        for row in rows:
            # print(row)
            primary_accession = row['primary_accession']
            host = row['host']
            alignment = row['alignment']
            product = row['product']
            cds_start = row['cds_start']
            cds_end = row['cds_end']

            if alignment is None:
                continue
            sub_seq = alignment[int(cds_start): int(cds_end)+1]
            cds_seq = sub_seq.replace("-", "")
            if len(cds_seq) < 3:
                continue
            protein = str(Seq(cds_seq).translate(to_stop=False))
            
            translated_results.append({
                "primary_accession": primary_accession,
                "host": host,
                "protein": protein,
                "region":product
            })

        results = {"reference_protein": translated_ref_results, "translated_sequences": translated_results}
        # print(translated_results)
        return results
    
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



    def get_mutations(self, codons, region, include_metadata=True, sequence_ids=None, hosts=None):
        """
        Fetch mutations for specific codons in a given genomic region across host samples.

        Args:
            hosts (list): List of host species to filter sequences.
            codons (list): List of codon indices to analyze.
            region (str): Genomic product/region name to restrict mutation analysis.

        Returns:
            list: A list of sequence alignment records with annotated codon-level mutations.
        
        Raises:
            ValueError: If required parameters (hosts, codons, or region) are missing.
        """
        if not region:
            raise ValueError("No region chosen")
        if not codons:
            raise ValueError("No codons chosen")

        self.sequence_ids = sequence_ids
        self.hosts = hosts
        self.codons = codons
        self.region = region
        alignments = self.__get_alignments()








        mutations = self.__parse_mutations(alignments)
        return mutations




    def get_mutation_regions_and_codons(self):
        """
        Get start and end positions of coding regions, annotated with codon labels.

        Returns:
            list: A list of feature dictionaries, each annotated with codon_start and codon_end.
        """
        master_reference = self.__get_master_reference()
        for feature in master_reference:
            codon_start, codon_end = get_codon_labeling(feature["cds_start"], feature["cds_end"])
            feature["codon_start"] = codon_start
            feature["codon_end"] = codon_end
            feature["proudct"] = feature["product"]
        return master_reference

    def __parse_mutations(self, alignments):
        """
        Annotate mutations by comparing aligned codons to the reference.

        Args:
            alignments (list): Aligned sequences from host species.
            master_reference (list): Genomic feature reference list.

        Returns:
            list: Annotated alignment records with mutations.
        """

        for alignment in alignments:
            ref_start = alignment["cds_start"]
            ref_end = alignment["cds_end"]
            alignment["mutations"] = {}
            sub_seq = alignment["alignment"][ref_start:ref_end+1]
            codons_new = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]

            for codon in self.codons:
                if 0 <= int(codon) < len(codons_new):
                    selected_codon = codons_new[int(codon)]
                    if "-" not in selected_codon:
                        if len(selected_codon) % 3 == 0:
                            protein = translateCodon(selected_codon)
                            alignment["mutations"][str(codon)] = protein
                    elif selected_codon != "---":
                        alignment["mutations"][str(codon)] = "X"

                else:
                    alignment["mutations"][str(codon)] = ""

        return alignments

    def __get_alignments(self):
        """
        Fetch sequence alignments for the given hosts.

        Returns:
            list: List of alignment and metadata records.
        """
        

        #If using sequence_ids -- do nothing 

        #if using other filters 

        formatted_hosts = ', '.join(['%s'] * len(self.hosts))
        query = f'''
            SELECT s.*, m.*, f.* 
            FROM sequence_alignment s 
            LEFT JOIN meta_data m ON s.sequence_id = m.primary_accession 
            LEFT JOIN features f ON s.sequence_id = f.accession 
            WHERE m.host IN ({formatted_hosts}) 
            AND f.product = '{self.region}';
        '''

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, self.hosts)
            alignments = dictfetchall(cursor)
        print(alignments)
        return alignments

    def __get_alignments2(self):
        """
        Fetch sequence alignments for the given hosts.

        Returns:
            list: List of alignment and metadata records.
        """
        formatted_hosts = ', '.join(['%s'] * len(self.hosts))
        query = f'''
            SELECT s.*, m.*, f.* 
            FROM sequence_alignment s 
            LEFT JOIN meta_data m ON s.sequence_id = m.primary_accession 
            LEFT JOIN features f ON s.sequence_id = f.accession 
            WHERE m.host IN ({formatted_hosts}) 
            AND f.product = '{self.region}';
        '''

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, self.hosts)
            alignments = dictfetchall(cursor)

        return alignments

    def __get_master_reference(self):
        """
        Retrieve reference feature annotations for the specified region.

        Returns:
            list: List of feature dictionaries.
        """
        with connections[self.database].cursor() as cursor:
            if hasattr(self, 'region'):
                cursor.execute(
                    "SELECT product, cds_start, cds_end FROM features WHERE accession=%s AND product=%s",
                    [self.reference_sequence, self.region]
                )
            else:
                cursor.execute(
                    "SELECT product, cds_start, cds_end FROM features WHERE accession=%s",
                    [self.reference_sequence]
                )
            master_reference = dictfetchall(cursor)

        return master_reference


