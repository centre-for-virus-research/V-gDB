from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

    
# from models import sequences_helpers as

class Taxonomy:

    def __init__(self, database, filters=None):

        self.database = database  
        self.filters = filters

    def _add_filters(self, params):

        where_clauses = []
        where_params = []
        for key, value in params.items():
            if isinstance(value, list):
                placeholders = ', '.join(['%s'] * len(value))
                where_clauses.append(f"{key} IN ({placeholders})")
                where_params.extend(value)
            else:
                where_clauses.append(f"{key} = %s")
                where_params.append(value)
        print(where_clauses, where_params)

        where_str = ' AND '.join(where_clauses)
        return where_str, where_params


    def get_phylum(self, params):
        where_clauses = ['phylum IS NOT NULL']
        filter_params = []
        if params:
            where_str, filter_params = self._add_filters(params)
            where_clauses.append(where_str)

        # query = 'SELECT DISTINCT(class) FROM host_lineage WHERE  ORDER BY class ASC'
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""
        
        query = f"""
            SELECT DISTINCT(phylum) FROM host_lineage
            {where_sql}
            ORDER BY phylum ASC
        """
        print(query, filter_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, filter_params)
            results = dictfetchall(cursor)
        return results

    def get_class(self, params):
        where_clauses = ['class IS NOT NULL']
        filter_params = []
        if params:
            where_str, filter_params = self._add_filters(params)
            where_clauses.append(where_str)

        # query = 'SELECT DISTINCT(class) FROM host_lineage WHERE  ORDER BY class ASC'
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        query = f"""
            SELECT DISTINCT(class) FROM host_lineage
            {where_sql}
            ORDER BY class ASC
        """
        print(query, filter_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, filter_params)
            results = dictfetchall(cursor)
        return results

    def get_order(self, params):
        where_clauses = ['order_category IS NOT NULL']
        filter_params = []
        if params:
            where_str, filter_params = self._add_filters(params)
            where_clauses.append(where_str)

        # query = 'SELECT DISTINCT(class) FROM host_lineage WHERE  ORDER BY class ASC'
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        query = f"""
            SELECT DISTINCT(order_category) FROM host_lineage
            {where_sql}
            ORDER BY order_category ASC
        """
        print(query, filter_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, filter_params)
            results = dictfetchall(cursor)
        return results

    def get_family(self, params):
        where_clauses = ['family IS NOT NULL']
        filter_params = []
        if params:
            where_str, filter_params = self._add_filters(params)
            where_clauses.append(where_str)

        # query = 'SELECT DISTINCT(class) FROM host_lineage WHERE  ORDER BY class ASC'
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        query = f"""
            SELECT DISTINCT(family) FROM host_lineage
            {where_sql}
            ORDER BY family ASC
        """
        print(query, filter_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, filter_params)
            results = dictfetchall(cursor)
        return results


    def get_genus(self, params):
        where_clauses = ['genus IS NOT NULL']
        filter_params = []
        if params:
            where_str, filter_params = self._add_filters(params)
            where_clauses.append(where_str)

        # query = 'SELECT DISTINCT(class) FROM host_lineage WHERE  ORDER BY class ASC'
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        query = f"""
            SELECT DISTINCT(genus) FROM host_lineage
            {where_sql}
            ORDER BY genus ASC
        """
        print(query, filter_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, filter_params)
            results = dictfetchall(cursor)
        return results

    def get_species(self, params):
        where_clauses = ['species IS NOT NULL']
        filter_params = []
        if params:
            where_str, filter_params = self._add_filters(params)
            where_clauses.append(where_str)

        # query = 'SELECT DISTINCT(class) FROM host_lineage WHERE  ORDER BY class ASC'
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        query = f"""
            SELECT DISTINCT(species), taxa_id FROM host_lineage
            {where_sql}
            ORDER BY species ASC
        """
        print(query, filter_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, filter_params)
            results = dictfetchall(cursor)
            taxa_ids = [item["taxa_id"] for item in results]
            placeholders = ', '.join(['%s'] * len(taxa_ids))
            common_taxa = f"taxa_id IN ({placeholders})"
            query2 = f""" SELECT DISTINCT(name) AS species 
                            FROM host_taxa 
                            WHERE {common_taxa}
                        """

            cursor.execute(query2, taxa_ids)
            results2 = dictfetchall(cursor)

            combined_results = results + results2
            
        return combined_results