from django.db import connections
import csv
from models.helpers import *

class Sequences:
    """
    A class to handle operations related to sequence metadata, alignments, 
    and filtered retrieval from a database in a Django project.
    """

    def __init__(self, database, filters=None):
        """
        Initialize the Sequences instance.

        Args:
            database (str): The name of the database alias to connect to.
            primary_accession (str, optional): Identifier for a specific sequence. Defaults to None.
            filters (dict, optional): Dictionary of filters for querying metadata. Defaults to None.
        """
        self.database = database  
        self.filters = filters

    def get_sequences_meta_data(self):
        """
        Retrieve all sequence metadata records from the database.

        Returns:
            list: A list of dictionaries containing metadata records.
        """
        with connections[self.database].cursor() as cursor:
            cursor.execute('SELECT * FROM meta_data ORDER BY create_date DESC;')
            result = dictfetchall(cursor)
        
        return result
    
    def get_sequence_meta_data(self, primary_accession):
        """
        Retrieve metadata and alignment details for a specific sequence, 
        including insertions, features, and reference sequence.

        Returns:
            dict: A dictionary containing metadata and alignment details.

        Raises:
            ValueError: If `primary_accession` is not provided or not found in the database.
        """
        result = {}

        if not primary_accession:
            raise ValueError("Primary accession can not be blank")

        with connections[self.database].cursor() as cursor:

            # Get main metadata
            cursor.execute("SELECT * FROM meta_data WHERE primary_accession = %s;", [primary_accession])
            sequence = dictfetchall(cursor)
            
            if not sequence:
                raise ValueError(f"Sequence with primary_accession {primary_accession} not found")
            
            result["meta_data"] = sequence[0]

            # Get alignment information
            cursor.execute("SELECT * FROM sequence_alignment WHERE sequence_id=%s", [primary_accession])
            alignment = dictfetchall(cursor)

            if alignment:
                result["alignment"] = alignment[0]

                # Add insertions
                cursor.execute("SELECT * FROM insertions WHERE accession = %s", [primary_accession])
                result["alignment"]["insertions"] = dictfetchall(cursor)

                # Add features (hardcoded accession in original code)
                cursor.execute("SELECT * FROM features WHERE accession = %s ORDER BY cds_start", ['EF437215'])
                result["alignment"]["features"] = dictfetchall(cursor)

                # Get reference sequence
                cursor.execute("SELECT sequence FROM sequences WHERE header = %s", [result["alignment"]["alignment_name"]])
                result["alignment"]["ref_seq"] = dictfetchall(cursor)[0]["sequence"]

            # Get regional info if country exists
            if result["meta_data"].get("country"):
                cursor.execute(
                    "SELECT * FROM m49_country WHERE full_name=%s",
                    [result["meta_data"]["country"].split(':')[0]]
                )
                country = dictfetchall(cursor)
                if country:
                    result["meta_data"]["region"] = country[0]

        return result

    def get_sequences_meta_data_by_filters(self):
        """
        Retrieve metadata based on filters defined in the `filters` dictionary.

        Returns:
            list: A list of dictionaries matching the filter criteria.
        """
        if not self.filters:
            # No filters provided, return all metadata
            return self.get_sequences_meta_data()

        where_clauses = []
        params = []

        def add_filter_clause(key, value, operator='='):
            """Add WHERE clause for numeric or comparison-based filters."""
            where_clauses.append(f"{key} {operator} %s")
            params.append(value)

        def add_filter_in_clause(key, value):
            """Add simple equality filter for strings/other fields."""
            where_clauses.append(f"{key} = %s")
            params.append(value)

        for key, value in self.filters.items():
            if key == 'length_lower':
                add_filter_clause('length', value, operator='>=')
            elif key == 'length_upper':
                add_filter_clause('length', value, operator='<=')
            else:
                add_filter_in_clause(key, value)

        where_str = ' AND '.join(where_clauses)
        query = f"SELECT * FROM meta_data WHERE {where_str} ORDER BY create_date DESC;"

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, params)
            result = dictfetchall(cursor)

        return result

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
            codons = get_codon_labeling(feature["cds_start"], feature["cds_end"])
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

