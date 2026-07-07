from .helpers import *
from django.db import connections

def _get_aligned_sequences_and_features_from_reference(cursor, reference_accession):

    query = f"""
                SELECT * 
                FROM sequence_alignment sa
                WHERE sa.alignment_name = %s
            """
                  # JOIN features f ON f.accession = sa.alignment_name
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
    data = fetch_all(cursor, query, params=None)
    return data

def _get_m49_country_data(cursor):
    query = "SELECT display_name, id, m49_code FROM m49_country"
    data = fetch_all(cursor, query, params=None)
    return data

def _get_taxa_from_host_taxa_id(cursor, host_taxa_id):
    query = "SELECT * FROM host_lineage WHERE taxa_id=%s"
    params = [host_taxa_id]
    data = fetch_one(cursor, query, params)
    return data

def _add_standard_filters(key, value, exclude):

    params, where_clauses, placeholders = [], [], []
    comparison = get_comparison(value, exclude)

    if isinstance(value, list):
        placeholders = ', '.join(['%s'] * len(value))
        where_clauses = f"{key} {comparison} ({placeholders})"
        params.extend(value)
    else:
        where_clauses = f"{key} {comparison} %s"
        params.append(value)

    return where_clauses, params

def _add_exclude_clade_filters(key, major, minor=None, exclude=True):

    params, where_clauses, placeholders_major, placeholders_minor = [], [], [], []
    # comparison = get_comparison(value, exclude)

    if isinstance(major, list):
        placeholders_major = ', '.join(['%s'] * len(major))
        where_clauses_major = f"EPA_major_clade IN ({placeholders_major})"
        params.extend(major)
    else:
        where_clauses_major = f"EPA_major_clade = %s"
        params.append(major)

    if minor:
        if isinstance(minor, list):
            placeholders_minor = ', '.join(['%s'] * len(minor))
            where_clauses_minor = f"EPA_minor_clade IN ({placeholders_minor})"
            params.extend(minor)
        else:
            where_clauses_minor = f"EPA_minor_clade = %s"
            params.append(minor)

    if minor:
        where_clauses = f"NOT ({where_clauses_major} AND {where_clauses_minor})"
    else:
        where_clauses = f"NOT ({where_clauses_major})"
        

    return where_clauses, params

def _add_genome_coverage_filters():
    return

def _add_region_filters(clauses, comparison):

    region_where_str = ' AND '.join(clauses)
    region_sql = f"""
                country_validated {comparison} (
                    SELECT m49_code
                    FROM m49_country 
                    WHERE {region_where_str})
                """
    return region_sql

    # where_clauses, params = _add_standard_filters("display_name", value, exclude)

    # region_sql = f"""
    #             country_validated IN (
    #                 SELECT m49_code
    #                 FROM m49_country 
    #                 WHERE {where_clauses})
    #             """
    
    # return region_sql, params
def _add_taxonomy_host_filters(value, comparison):
    # taxa_where_str = ' AND '.join(clauses)
    params = []
    common_clause, common_param = _add_standard_filters("common_name", value, False)

    common_sql = f""" SELECT taxa_id FROM host_taxa WHERE {common_clause} """

    sql = f""" ( host_taxa_id {comparison} ( {common_sql} ) 
                )
            """
    
    params.extend(common_param)
    print("COMMON ", common_sql, common_param)
    print("-----")
    recursive_sql = f""" WITH RECURSIVE taxa_tree(id) AS (
                    SELECT parent_taxa_id
                    FROM host_children
                    WHERE parent_taxa_id IN ({common_sql})

                UNION

                    SELECT hc.child_taxa_id
                    FROM host_children hc
                    JOIN taxa_tree tt
                    ON hc.parent_taxa_id = tt.id 
                )
            SELECT DISTINCT id FROM taxa_tree;
            """
    print(recursive_sql)
    with connections["RABV"].cursor() as cursor:
        cursor.execute(recursive_sql, params)
        rows = cursor.fetchall()
        ids = [row[0] for row in rows]

    print("IDS", ids)

    # return recursive_sql, params
    if len(ids) == 0:
        clause = sql 
        param = params
    else:
        clause, param = _add_standard_filters("host_taxa_id", ids, False)

    return clause, param

def _add_taxonomy_species_filters(value, comparison):
    # taxa_where_str = ' AND '.join(clauses)
    params = []
    common_clause, common_param = _add_standard_filters("host", value, False)
    scientific_clause, scientific_param = _add_standard_filters("species", value, False)

    common_sql = f""" SELECT host_taxa_id FROM meta_data WHERE {common_clause} """
    scientific_sql = f""" SELECT taxa_id FROM host_lineage WHERE {scientific_clause} """

    sql = f""" ( host_taxa_id {comparison} ( {common_sql} ) 
                OR host_taxa_id {comparison} ({scientific_sql})
                )
            """
    
    params.extend(common_param)
    params.extend(scientific_param)

    recursive_sql = f""" WITH RECURSIVE taxa_tree(id) AS (
                    SELECT parent_taxa_id
                    FROM host_children
                    WHERE parent_taxa_id IN ({common_sql}) OR parent_taxa_id IN ({scientific_sql})

                UNION

                    SELECT hc.child_taxa_id
                    FROM host_children hc
                    JOIN taxa_tree tt
                    ON hc.parent_taxa_id = tt.id 
                )
            SELECT DISTINCT id FROM taxa_tree;
            """
    print(recursive_sql)
    with connections["RABV"].cursor() as cursor:
        cursor.execute(recursive_sql, params)
        rows = cursor.fetchall()
        ids = [row[0] for row in rows]

    print("IDS", ids)

    # return recursive_sql, params
    if len(ids) == 0:
        clause = sql 
        param = params
    else:
        clause, param = _add_standard_filters("host_taxa_id", ids, False)

    return clause, param


def _add_taxonomy_filters(clauses, comparison):
    taxa_where_str = ' AND '.join(clauses)
    taxa_sql = f""" host_taxa_id {comparison} (
                SELECT taxa_id
                FROM host_lineage
                WHERE {taxa_where_str}
            )
        """
    return taxa_sql

def get_comparison(value, exclude):
    if isinstance(value, list):
        comparison = "IN"
        if (exclude):
            comparison = "NOT IN"
    else:
        comparison = "="
        if (exclude):
            comparison = "!="
    return comparison


def recursive_taxa_search(ids):

    sql = f""" WITH RECURSIVE taxa_tree(id) AS (
                    SELECT parent_taxa_id
                    FROM host_children
                    WHERE parent_taxa_id IN (9822)

                UNION

                    SELECT hc.child_taxa_id
                    FROM host_children hc
                    JOIN taxa_tree tt
                    ON hc.parent_taxa_id = tt.id 
                )
            SELECT DISTINCT id FROM taxa_tree;
            """