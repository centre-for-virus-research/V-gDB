from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

    
from models import sequences_helpers as sh


def _get_tree_from_type(cursor, tree_type, segment):

    if segment:
        query = "SELECT * FROM tree WHERE tree_type = %s AND segment = %s;"
        params = [tree_type, segment]
    else:
        query = "SELECT * FROM tree"
        params = None
    
    data = fetch_one(cursor, query, params)
    
    return data

def _get_meta_data_for_tree(cursor, segment):
    if segment:
        query = f"""SELECT Parsed_strain, host, serotype
                    FROM meta_data WHERE segment = %s;""" 

        
        params = [segment]
    else:
        query = f"""SELECT md.primary_accession, md.EPA_major_clade, md.EPA_minor_clade, md.collection_year, c.display_name  as country, h.*
                    FROM meta_data md 
                    JOIN m49_country c ON c.m49_code = md.country_validated
                    JOIN host_lineage h ON h.taxa_id = md.host_taxa_id
                """
        params = None
    data = fetch_all(cursor, query, params)

    return data

class Phylogeny:

    def __init__(self, database, filters=None):

        self.database = database  
        self.filters = filters

    def get_tree(self):

        tree_type=None
        segment=None

        if self.filters:
            tree_type = self.filters["tree_type"]
            segment = self.filters["segment"]

        results = {}
        with connections[self.database].cursor() as cursor:

            tree = _get_tree_from_type(cursor, tree_type, segment)
            meta_data = _get_meta_data_for_tree(cursor, segment)

            
            results = {
                        "tree": tree,
                        "meta_data": meta_data
                        }   

        return results
