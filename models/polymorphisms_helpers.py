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

def _get_polymorphim_sequence(cursor, id):

    query = f"""
                SELECT *
                FROM completed_signatures_only
                WHERE primary_accession = %s;
            """
    params = [id]
    results = fetch_all(cursor, query, params)

    return results

def _get_mutation_prevalence(cursor, id):

    accession_query = f""" 
                            SELECT cso.primary_accession
                            FROM completed_signatures_only cso
                            WHERE cso.signature_id = %s
                        """

    reference_query = f""" 
                SELECT accession, reference_accession, cds_start, cds_end
                FROM features
                WHERE product LIKE %s
                AND accession IN ({accession_query})
            """
    
#     reference_query = f"""
# SELECT f.accession, f.reference_accession, f.cds_start, f.cds_end
# FROM features f
# WHERE f.product LIKE "%NS3%"
#                 AND f.accession IN ( 
#                         SELECT cso.primary_accession
#                         FROM completed_signatures_only cso
#                         WHERE cso.signature_id = "NS3:117L"
#                     )
#                         """
    # meta_data_query = f"""
    #         SELECT sa.header as sequence_id, sa.sequence as alignment, f.reference_accession, f.cds_start, f.cds_end, md.host, md.nearest_reference_genotype, md.nearest_reference_subtype, md.host_scientific_name, md.country
    #         FROM sequences sa 
    #         LEFT JOIN features f ON f.accession = sa.header
    #         LEFT JOIN meta_data md ON md.primary_accession = sa.header
    #         WHERE f.product LIKE %s
    #         AND sa.header IN ({accession_query})
    #     """
    meta_data_query = f"""
                    SELECT sa.sequence_id, sa.alignment, f.reference_accession, f.cds_start, f.cds_end, md.host, md.nearest_reference_genotype, md.nearest_reference_subtype, md.host_scientific_name, md.country
                    FROM sequence_alignment sa 
                    LEFT JOIN features f ON f.accession = sa.sequence_id
                    LEFT JOIN meta_data md ON md.primary_accession = sa.sequence_id
                    WHERE f.product LIKE %s
                    AND sa.sequence_id IN ({accession_query})
                """
    
    clades_count_query = f"""
                        SELECT nearest_reference_genotype, COUNT(nearest_reference_genotype) 
                        FROM meta_data
                        WHERE exclusion_status = 0
                        GROUP BY nearest_reference_genotype
                    """
    
    genotype_count_query = f"""
                        SELECT nearest_reference_genotype, COUNT(nearest_reference_genotype) as count
                        FROM meta_data
                        WHERE exclusion_status = 0
                        GROUP BY nearest_reference_genotype
                    """
    subtype_count_query = f"""
                        SELECT nearest_reference_genotype, nearest_reference_subtype, COUNT(nearest_reference_subtype) as count
                        FROM meta_data
                        WHERE exclusion_status = 0
                        GROUP BY nearest_reference_genotype, nearest_reference_subtype
                    """
       
    params_reference = [f"%NS3%", id]
    params_md = [f"%NS3%", id]
    results_reference = fetch_all(cursor, reference_query, params_reference)
    results_meta_data = fetch_all(cursor, meta_data_query, params_md)
    results_genotype_count = fetch_all(cursor, genotype_count_query, params=None)
    results_subtype_count = fetch_all(cursor, subtype_count_query, params=None)
    # results_meta_data = None

    return {"reference":results_reference, 
            "meta_data":results_meta_data, 
            "genotype_count": results_genotype_count, 
            "subtype_count": results_subtype_count}