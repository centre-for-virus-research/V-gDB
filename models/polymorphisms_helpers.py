from .helpers import *
from django.db import connections

def _get_all_polymorphims(cursor, filters=None):
    where_clauses = []
    where_str = []
    params = None
    if filters:
        params = []
        # where_str, filter_params = self._add_filters()
        print(filters.items())
        for key, value in filters.items():
            value = value.split(',')
            print("VALUE", value)
            if key.startswith("exclude_"): 
                key=key[8:]
            print("KEY ", key, "VALUE ", value)
            if isinstance(value, list):
                comparison = "IN"
                if ("exclude_"+key in filters):
                    comparison = "NOT IN"

                placeholders = ', '.join(['%s'] * len(value))
                where_clauses.append(f"{key} {comparison} ({placeholders})")
                params.extend(value)
            else:
                comparison = "="
                if ("exclude_"+key in filters):
                    comparison = "!="
                where_clauses.append(f"{key} {comparison} %s")
                params.append(value)

            # where_clauses.append(where_str)
    print("WHERE clauses", where_clauses)
    filter_where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""


    query = f"""
                SELECT signature_id, signature_kind, protein_name, aa_position, mutation_type, drug
                FROM mutation_catalog
                {filter_where_sql};
            """
    results = fetch_all(cursor, query, params)

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
    print(id, "ID")
    gene = id.split(":")[0]
    print("GENE", gene)

    accession_query = f""" 
                            SELECT cso.primary_accession
                            FROM completed_signatures_only cso
                            WHERE cso.signature_id = %s
                        """

    reference_query = f"""
                        SELECT accession, cds_start, cds_end
                        FROM features
                        WHERE product LIKE %s
                        AND accession IN (
                            SELECT DISTINCT reference_accession
                            FROM features
                            WHERE accession IN  ({accession_query}))
                    """

    meta_data_query = f"""
                    SELECT sa.sequence_id, sa.alignment, f.reference_accession, md.host, md.EPA_major_clade, md.EPA_minor_clade, md.host_scientific_name, md.country
                    FROM sequence_alignment sa 
                    LEFT JOIN features f ON f.accession = sa.sequence_id
                    LEFT JOIN meta_data md ON md.primary_accession = sa.sequence_id
                    WHERE f.product LIKE %s
                    AND sa.sequence_id IN ({accession_query})
                """
       
    genotype_count_query = f"""
                        SELECT EPA_major_clade, COUNT(EPA_major_clade) as count
                        FROM meta_data
                        WHERE exclusion_status = 0
                        GROUP BY EPA_major_clade
                    """
    subtype_count_query = f"""
                        SELECT EPA_major_clade, EPA_minor_clade, COUNT(EPA_minor_clade) as count
                        FROM meta_data
                        WHERE exclusion_status = 0
                        GROUP BY EPA_major_clade, EPA_minor_clade
                    """
    params_reference = [f"%{gene}%", id]
    params_md = [f"%{gene}%", id]
    results_reference = fetch_all(cursor, reference_query, params_reference)
    results_meta_data = fetch_all(cursor, meta_data_query, params_md)
    results_genotype_count = fetch_all(cursor, genotype_count_query, params=None)
    results_subtype_count = fetch_all(cursor, subtype_count_query, params=None)
    # results_meta_data = None

    print(results_reference)

    return {

            "reference":results_reference, 
            "meta_data":results_meta_data,
            "genotype_count": results_genotype_count, 
            "subtype_count": results_subtype_count

            }

# def _get_mutation_prevalence(cursor, id):

#     accession_query = f""" 
#                             SELECT cso.primary_accession
#                             FROM completed_signatures_only cso
#                             WHERE cso.signature_id = %s
#                         """

#     # reference_query = f""" 
#     #             SELECT accession, reference_accession, cds_start, cds_end
#     #             FROM features
#     #             WHERE product LIKE %s
#     #             AND accession IN ({accession_query})
#     #         """

#     reference_query = f"""
#                         SELECT accession, cds_start, cds_end
#                         FROM features
#                         WHERE product LIKE %s
#                         AND reference_accession IN (select reference_accession from features WHERE accession IN ({accession_query}) )

#     """
    
# #     reference_query = f"""
# # SELECT f.accession, f.reference_accession, f.cds_start, f.cds_end
# # FROM features f
# # WHERE f.product LIKE "%NS3%"
# #                 AND f.accession IN ( 
# #                         SELECT cso.primary_accession
# #                         FROM completed_signatures_only cso
# #                         WHERE cso.signature_id = "NS3:117L"
# #                     )
# #                         """
#     # meta_data_query = f"""
#     #         SELECT sa.header as sequence_id, sa.sequence as alignment, f.reference_accession, f.cds_start, f.cds_end, md.host, md.EPA_major_clade, md.EPA_minor_clade, md.host_scientific_name, md.country
#     #         FROM sequences sa 
#     #         LEFT JOIN features f ON f.accession = sa.header
#     #         LEFT JOIN meta_data md ON md.primary_accession = sa.header
#     #         WHERE f.product LIKE %s
#     #         AND sa.header IN ({accession_query})
#     #     """
#     meta_data_query = f"""
#                     SELECT sa.sequence_id, sa.alignment, f.reference_accession, f.cds_start, f.cds_end, md.host, md.EPA_major_clade, md.EPA_minor_clade, md.host_scientific_name, md.country
#                     FROM sequence_alignment sa 
#                     LEFT JOIN features f ON f.accession = sa.sequence_id
#                     LEFT JOIN meta_data md ON md.primary_accession = sa.sequence_id
#                     WHERE f.product LIKE %s
#                     AND sa.sequence_id IN ({accession_query})
#                 """
    
#     clades_count_query = f"""
#                         SELECT EPA_major_clade, COUNT(EPA_major_clade) 
#                         FROM meta_data
#                         WHERE exclusion_status = 0
#                         GROUP BY EPA_major_clade
#                     """
    
#     genotype_count_query = f"""
#                         SELECT EPA_major_clade, COUNT(EPA_major_clade) as count
#                         FROM meta_data
#                         WHERE exclusion_status = 0
#                         GROUP BY EPA_major_clade
#                     """
#     subtype_count_query = f"""
#                         SELECT EPA_major_clade, EPA_minor_clade, COUNT(EPA_minor_clade) as count
#                         FROM meta_data
#                         WHERE exclusion_status = 0
#                         GROUP BY EPA_major_clade, EPA_minor_clade
#                     """
       
#     params_reference = [f"%NS3%", id]
#     params_md = [f"%NS3%", id]
#     results_reference = fetch_all(cursor, reference_query, params_reference)
#     results_meta_data = fetch_all(cursor, meta_data_query, params_md)
#     results_genotype_count = fetch_all(cursor, genotype_count_query, params=None)
#     results_subtype_count = fetch_all(cursor, subtype_count_query, params=None)
#     # results_meta_data = None

#     print(results_reference)

#     return {"reference":results_reference, 
#             "meta_data":results_meta_data, 
#             "genotype_count": results_genotype_count, 
#             "subtype_count": results_subtype_count}