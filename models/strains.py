from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

    
from models import sequences_helpers as sh
from models import polymorphisms_helpers as ph

from models.filter_helper import FilterHelper

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

    def __build_pagination_query(self, where_clauses, filter_params, next_cursor, prev_cursor, items_per_page):

        pagination_clauses = where_clauses.copy()
        pagination_params = filter_params.copy()

        order_by = "ORDER BY strain_accession"
        limit = "LIMIT %s"
        pagination_params.append(items_per_page)

        if next_cursor:
            pagination_clauses.append("strain_accession > %s")
            pagination_params.insert(-1, next_cursor)

        elif prev_cursor:
            order_by = "ORDER BY strain_accession DESC"
            if prev_cursor != "0":
                pagination_clauses.append("strain_accession < %s")
                pagination_params.insert(-1, prev_cursor)

        pagination_where_sql = (
            f"WHERE {' AND '.join(pagination_clauses)}"
            if pagination_clauses else ""
        )

        query = f"""
                    SELECT *
                    FROM strain_segment_linkage
                    {pagination_where_sql}
                    {order_by}
                    {limit}
                """
        
        return query, pagination_params


    def get_strains(self, next_cursor, prev_cursor, items_per_page):
        
        where_clauses, filter_params = [], []

        print('DATABASE', self.database) #TODO: remove this once testing is finished


        ############# BUILD WHERE string if there are filters
        if self.filters:
            filterHelper = FilterHelper(filters = self.filters, database=self.database)
            where_str, filter_params = filterHelper.add_filters()
            # where_str, filter_params = self._add_filters_refactored()
            where_clauses.append(where_str)
        filter_where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""
        #############


        ############ BUILDING THE PAGINATION QUERY ############
        query, pagination_params = self.__build_pagination_query(where_clauses, filter_params, next_cursor, prev_cursor, items_per_page)
        print("PAGINATION QUERY", query, pagination_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, pagination_params)
            results = dictfetchall(cursor)

        if prev_cursor:
            results.reverse()

        total_count = self.__get_strains_count(filter_where_sql, filter_params)

        next_cursor_new, prev_cursor_new = None, None
        if results: 
            next_cursor_new = results[-1]["strain_accession"]
            prev_cursor_new = results[0]["strain_accession"]

        results_dict = {
                            "data": results,
                            "total_count": total_count,
                            "next_cursor": next_cursor_new,
                            "prev_cursor": prev_cursor_new,
                        }
        
        return results_dict


        
    def get_strain(self, parsed_strain):
        """
        Retrieve metadata and alignment details for a specific sequence, 
        including insertions, features, and reference sequence.

        Returns:
            dict: A dictionary containing metadata and alignment details.

        Raises:
            ValueError: If `primary_accession` is not provided or not found in the database.
        """
        result = {}

        print(parsed_strain)

        if not parsed_strain:
            raise ValueError("Strain can not be blank")

        with connections[self.database].cursor() as cursor:
            query = "SELECT * from meta_data where parsed_strain = %s"

            strain = fetch_all(cursor, query, [parsed_strain])
            if strain: 
                result["strain"] = strain

            
            # meta_data = sh._get_meta_data_from_primary_accession(cursor, primary_accession)
            
            # if not meta_data:
            #     raise ValueError(f"Sequence with primary_accession {primary_accession} not found")
            
            # result["meta_data"] = meta_data

            # sequence = sh._get_sequence_from_primary_accession(cursor, primary_accession)
            # result.update(sequence)

            # # Get regional info if country exists
            # if meta_data["country_validated"]:
            #     country_code = int(meta_data["country_validated"])
            #     regions = sh._get_region_from_country_code(cursor, country_code)
            #     result["regions"] = regions

            # # Add hose data
            # if meta_data["host_taxa_id"]:
            #     host_taxa_id = meta_data["host_taxa_id"]
            #     if self.database == 'FLU' or self.database == 'HCV':
            #         taxanomic_info = None
            #     else:
            #         taxanomic_info = sh._get_taxa_from_host_taxa_id(cursor, host_taxa_id)

            #     if taxanomic_info:
            #         result["taxanomic_info"] = taxanomic_info

            # query_alignment_dict = sh._get_query_alignment_from_primary_accession(cursor, primary_accession)
        
            # if query_alignment_dict:

            #     reference_accession = query_alignment_dict["alignment_name"]
                
            #     query_alignment_sequence = query_alignment_dict["alignment"].upper()
            #     if self.database == 'HCV' or self.database == 'FLU':

            #         insertions = sh._get_insertions_from_primary_accession(
            #             cursor,
            #             primary_accession
            #         )

            #         insertions = insertions[0]["insertion"] if len(insertions) > 0 else None

            #         # result["insertions"] = insertions if insertions else None
            #     else:
            #         insertions = query_alignment_dict["insertion"]
            #     result["insertions"] = insertions if insertions else None

            #     # Get aligned reference sequence
            #     reference_alignment_dict = sh._get_query_alignment_from_primary_accession(cursor, reference_accession)
            #     reference_alignment_sequence = reference_alignment_dict["alignment"]

                
            #     # Get features
            #     features = sh._get_features_from_primary_accession(cursor, primary_accession)

            #     alignment = {
            #                     "reference_accession": reference_accession,
            #                     "query_alignment_sequence":query_alignment_sequence,
            #                     "reference_alignment_sequence":reference_alignment_sequence,
            #                     "features": features,
            #                 }    
            
            #     result["alignment"] = alignment

            # if self.database == 'HCV':
            #     mutations = ph._get_polymorphim_sequence(cursor, primary_accession)
            #     if mutations:
            #         result["mutations"] = mutations

        return result


    def __get_strains_count(self, where_sql, params):

        query = f"""
                    SELECT COUNT(strain_accession)
                    FROM strain_segment_linkage
                    {where_sql}
                """

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, params)
            total_count = cursor.fetchone()[0]

        return total_count