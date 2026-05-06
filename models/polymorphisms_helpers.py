from .helpers import *
from django.db import connections

def _get_all_polymorphims(cursor):

    query = f"""
                SELECT signature_id, COUNT(*) as signature_count, signature_kind
                FROM mutation_catalog 
                GROUP BY signature_id;
            """
    results = fetch_all(cursor, query, params=None)

    return results

def _get_polymorphim(cursor, id):

    query = f"""
                SELECT *
                FROM mutation_catalog 
                WHERE signature_id=%s;
            """
    params = [id]
    results = fetch_all(cursor, query, params)

    return results

def _get_polymorphim_sequences(cursor, id):

    query = f"""
                SELECT *
                FROM completed_signatures_only
                WHERE signature_id = %s;
            """
    params = [id]
    results = fetch_all(cursor, query, params)

    return results
