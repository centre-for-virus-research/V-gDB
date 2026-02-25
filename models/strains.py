from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

from models import strains_helpers as sh
from models import sequences_helpers as seqh

class Strains:
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

    def get_strains(self):

        
        with connections[self.database].cursor() as cursor:
            result = sh._get_all_strains(cursor)

        return result

    def get_strain(self, strain_id):

        if not strain_id:
            raise ValueError("Strain can not be blank")

        segment_data = {}

        with connections[self.database].cursor() as cursor:

            segment_data = sh._get_segment_meta_data_from_strain_id(cursor, strain_id)

            if segment_data:
                for segment in segment_data:
                    primary_accession = segment["primary_accession"]

                    insertions = seqh._get_insertions_from_primary_accession(cursor, primary_accession)
                    if insertions: 
                        segment["insertions"] = insertions

                    query_alignment_dict = seqh._get_query_alignment_from_primary_accession(cursor, primary_accession)

                    if query_alignment_dict:
                        reference_accession = query_alignment_dict["alignment_name"]
                        query_alignment_sequence = query_alignment_dict["alignment"]

                        
                    
                        features = seqh._get_features_from_primary_accession(cursor, primary_accession)
                        reference_alignment_dict = seqh._get_query_alignment_from_primary_accession(cursor, reference_accession)
                        reference_alignment_sequence = reference_alignment_dict["alignment"]

                        segment["query_alignment_sequence"] = query_alignment_sequence
                        segment["reference_alignment_sequence"] = reference_alignment_sequence
                        segment["features"] = features
            

        
        return segment_data
