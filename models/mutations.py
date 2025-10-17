from django.db import connections
from models.helpers import * 


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


