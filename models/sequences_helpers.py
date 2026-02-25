from .helpers import *

def _get_aligned_sequences_and_features_from_reference(cursor, reference_accession):

    query = f"""
                SELECT * 
                FROM sequence_alignment sa
                JOIN features f ON f.accession = sa.alignment_name
                WHERE sa.alignment_name = %s
            """
    params = [reference_accession]
    results = fetch_all(cursor, query, params)

    return results

def _get_meta_data_from_primary_accession(cursor, primary_accession):

    query = "SELECT * FROM meta_data WHERE primary_accession = %s;"
    params = [primary_accession]
    meta_data = fetch_one(cursor, query, params)
    
    return meta_data

def _get_sequence_from_primary_accession(cursor, primary_accession):

    query = "SELECT sequence FROM sequences WHERE header=%s"
    params = [primary_accession]

    sequence = fetch_one(cursor, query, params)
    return sequence

def _get_features_from_primary_accession(cursor, primary_accession):
    query = f"""
                SELECT * 
                FROM features 
                WHERE accession=%s 
                ORDER BY cds_start
            """
    params = [primary_accession]
    features = fetch_all(cursor, query, params)

    return features

def _get_insertions_from_primary_accession(cursor, primary_accession):
    query = "SELECT * FROM insertions WHERE accession = %s"
    params = [primary_accession]

    insertions = fetch_all(cursor, query, params)
    return insertions

def _get_query_alignment_from_primary_accession(cursor, primary_accession):

    query = "SELECT * FROM sequence_alignment WHERE sequence_id=%s"
    params = [primary_accession]

    query_alignment = fetch_one(cursor, query, params)
    return query_alignment

def _get_region_from_country_code(cursor, country_code):

    query = "SELECT * FROM m49_country WHERE m49_code=%s"
    params = [country_code]

    regions = fetch_one(cursor, query, params)
    return regions

def _get_country_meta_data(cursor):
    query = "SELECT country_validated FROM meta_data WHERE country_validated IS NOT NULL"
    data = fetchall(cursor, query, params=None)
    return data

def _get_m49_country_data(cursor):
    query = "SELECT display_name, id, m49_code FROM m49_country"
    data = fetchall(cursor, query, params=None)
    return data

def _get_taxa_from_host_taxa_id(cursor, host_taxa_id):
    query = "SELECT * FROM host_lineage WHERE taxa_id=%s"
    params = [host_taxa_id]
    data = fetch_one(cursor, query, params)
    return data