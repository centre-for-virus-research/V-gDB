from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

    
from models import polymorphisms_helpers as ph

class Polymorphisms:
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

    def get_polymorphisms(self, params):

        
        with connections[self.database].cursor() as cursor:
            result = ph._get_all_polymorphims(cursor, params)

        return result
    

    def get_polymorphism(self, id):

        
        with connections[self.database].cursor() as cursor:
            polymorphism = ph._get_polymorphim(cursor, id)

            sequences = ph._get_polymorphim_sequences(cursor, id)
            chart_data = ph._get_mutation_prevalence(cursor, id)


        result = {
                    "polymorphism": polymorphism,
                    "sequences": sequences,
                    "chart_data": chart_data
                }
        return result

    def get_polymorphisms_regions(self):

        
        with connections[self.database].cursor() as cursor:
            query = f"""
                SELECT *
                FROM features
                WHERE accession = 'NC_004102'
                AND product IN (%s, %s, %s);
            """
        
            features = fetch_all(cursor, query, ['nonstructural protein NS5A', 'protease/helicase protein NS3', 'RNA-dependent RNA polymerase NS5B'])

            query = f"""
                    SELECT sequence
                    FROM sequences
                    WHERE header = %s;
                """
            
            reference_sequence = fetch_one(cursor, query, ['NC_004102'])



        print(reference_sequence)

        result = {
                    "features": features,
                    "reference_sequence": reference_sequence['sequence']
                }
        return result
    
    

    