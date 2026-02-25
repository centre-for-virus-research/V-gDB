from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

    
from models import sequences_helpers as sh

class Sequences:
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

    def get_sequences(self, next_cursor, prev_cursor, items_per_page):
        where_clauses = []
        params = []
        filter_params = []

        if self.filters:
            where_str, filter_params = self._add_filters()
            where_clauses.append(where_str)

        filter_where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        count_query = f"""
            SELECT COUNT(primary_accession)
            FROM meta_data
            {filter_where_sql}
        """

        with connections[self.database].cursor() as cursor:
            cursor.execute(count_query, filter_params)
            total_count = cursor.fetchone()[0]

        pagination_clauses = where_clauses.copy()
        pagination_params = filter_params.copy()

        order_by = "ORDER BY primary_accession"
        limit = "LIMIT %s"
        pagination_params.append(items_per_page)

        if next_cursor:
            pagination_clauses.append("primary_accession > %s")
            pagination_params.insert(-1, next_cursor)

        elif prev_cursor:
            order_by = "ORDER BY primary_accession DESC"
            if prev_cursor != "0":
                pagination_clauses.append("primary_accession < %s")
                pagination_params.insert(-1, prev_cursor)

        pagination_where_sql = (
            f"WHERE {' AND '.join(pagination_clauses)}"
            if pagination_clauses else ""
        )

        query = f"""
            SELECT *
            FROM meta_data
            {pagination_where_sql}
            {order_by}
            {limit}
        """

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, pagination_params)
            results = dictfetchall(cursor)

        if prev_cursor:
            results.reverse()

        return {
            "data": results,
            "total_count": total_count,
            "next_cursor": results[-1]["primary_accession"] if results else None,
            "prev_cursor": results[0]["primary_accession"] if results else None,
        }

    def get_sequences_alignment(self, start_coordinate, end_coordinate, sequence_type):

        where_clauses = []
        params = []
        filter_params = []

        if self.filters:
            where_str, filter_params = self._add_filters()
            where_clauses.append(where_str)

        filter_where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        query = f"""
            SELECT primary_accession 
            FROM meta_data
            {filter_where_sql}
        """
        # print("QUERY: ", query)
        with connections[self.database].cursor() as cursor:

            # primary_accessions = fetch_all(cursor, query, filter_params)
            # accessions = [item['primary_accession'] for item in primary_accessions]

            # placeholders = ','.join(['?'] * len(accessions))

            query2 = f"""
                    SELECT sa.*, f.*
                    FROM sequence_alignment sa
                    JOIN features f ON f.accession = sa.sequence_id
                    WHERE sa.sequence_id IN ({query}) AND f.product = 'nucleoprotein N'
                    """
            # print(query2)
            cursor.execute(query2, filter_params)
            results = dictfetchall(cursor)

            alignments = self.__parse_alignments_new(results, start_coordinate, end_coordinate, sequence_type)
            
            # print(alignments)
            # rows = cursor.fetchall()

        

        return alignments
            # cursor.execute(query, filter_params)
            # data = fetchall()
            

    def __parse_alignments_new(self, alignments, start_coordinate, end_coordinate, sequence_type):
        results= []

        for a in alignments:
            if not start_coordinate and not end_coordinate:  # Use the full region coordinates

                ref_start = a["cds_start"]
                ref_end = a["cds_end"]

                if sequence_type == "codon":  # User wants a codon
                    codon_start, codon_end = get_codon_labeling(ref_start, ref_end)
            else:  # User chooses the coordinates themselves

                ref_start = int(start_coordinate)
                ref_end = int(end_coordinate)

                if sequence_type == "codon":  # User wants a codon
                    ref_start = a["cds_start"]
                    ref_end = a["cds_end"]
                    codon_start = int(start_coordinate)
                    codon_end = int(end_coordinate)

            sub_seq = a["alignment"][ref_start:ref_end+1]

            if (sequence_type == "codon"):
                codons = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]
                selected_codons = codons[codon_start-1:codon_end]
                sub_seq = ''.join(selected_codons)

            if set(sub_seq) != {"-"}:
                results.append(f">" + a["sequence_id"] + "\n" + sub_seq + "\n")

        return results
        

    def get_sequence(self, primary_accession):
        """
        Retrieve metadata and alignment details for a specific sequence, 
        including insertions, features, and reference sequence.

        Returns:
            dict: A dictionary containing metadata and alignment details.

        Raises:
            ValueError: If `primary_accession` is not provided or not found in the database.
        """
        result = {}

        if not primary_accession:
            raise ValueError("Primary accession can not be blank")

        with connections[self.database].cursor() as cursor:
            meta_data = sh._get_meta_data_from_primary_accession(cursor, primary_accession)
            
            if not meta_data:
                raise ValueError(f"Sequence with primary_accession {primary_accession} not found")
            
            result["meta_data"] = meta_data

            sequence = sh._get_sequence_from_primary_accession(cursor, primary_accession)
            result.update(sequence)

            # Get regional info if country exists
            if meta_data["country_validated"]:
                country_code = int(meta_data["country_validated"])
                regions = sh._get_region_from_country_code(cursor, country_code)
                result["regions"] = regions

            if meta_data["host_taxa_id"]:
                host_taxa_id = meta_data["host_taxa_id"]
                taxanomic_info = sh._get_taxa_from_host_taxa_id(cursor, host_taxa_id)

                if taxanomic_info:
                    result["taxanomic_info"] = taxanomic_info

            # Get insertions
            insertions = sh._get_insertions_from_primary_accession(cursor, primary_accession)
            if insertions:
                result["insertions"] = insertions

            query_alignment_dict = sh._get_query_alignment_from_primary_accession(cursor, primary_accession)
        
            if query_alignment_dict:

                reference_accession = query_alignment_dict["alignment_name"]
                query_alignment_sequence = query_alignment_dict["alignment"]

                # Get aligned reference sequence
                reference_alignment_dict = sh._get_query_alignment_from_primary_accession(cursor, reference_accession)
                reference_alignment_sequence = reference_alignment_dict["alignment"]
                
                # Get features
                features = sh._get_features_from_primary_accession(cursor, primary_accession)

                alignment = {
                                "reference_accession": reference_accession,
                                "query_alignment_sequence":query_alignment_sequence,
                                "reference_alignment_sequence":reference_alignment_sequence,
                                "features": features,
                            }    
            
                result["alignment"] = alignment

        return result

    def get_reference_sequence(self, primary_accession):

        result = {}

        # Grabbing reference features
        with connections[self.database].cursor() as cursor:

            features = sh._get_features_from_primary_accession(cursor, primary_accession)
        
        if not features:
            raise ValueError(f"Reference sequence with primary_accession {primary_accession} not found")
        
        
        for feature in features:
            codons = get_codon_labeling(feature["cds_start"], feature["cds_end"])
            feature["codon_start"] = codons[0]
            feature["codon_end"] = codons[1]
            
        result["features"] = features

        with connections[self.database].cursor() as cursor:

            aligned_sequences = sh._get_aligned_sequences_and_features_from_reference(cursor, primary_accession)
            filtered = [
                {k: d[k] for k in {"sequence_id", "alignment"} if k in d}
                for d in aligned_sequences
            ]

            key_map = {
                "sequence_id": "query_sequence_id",
                "alignment": "query_alignment_sequence"
            }

            filtered = [
                {new_key: d[old_key] for old_key, new_key in key_map.items() if old_key in d}
                for d in aligned_sequences
            ]

            reference_query = sh._get_query_alignment_from_primary_accession(cursor, primary_accession)
            reference_meta_data = sh._get_meta_data_from_primary_accession(cursor, primary_accession)
            
            result["query_aligned_sequences"] = filtered
            result["reference_alignment_sequence"] = reference_query["alignment"]
            result["reference_accession"] = primary_accession
            result["reference_meta_data"] = reference_meta_data

        return result

    def get_global_distribution_of_sequences(self):
        """
        Returns a list of unique countries (with m49 codes) and the number of sequences per country.
        """

        print(self.database)
        where_clauses = []
        params = []

        if self.filters:
            where_str, filter_params = self._add_filters()
            where_clauses.append(where_str)
            params.extend(filter_params)

        # Build query
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

                    
        query = f"""
                SELECT
                md.country_validated as m49_code,
                m.*,
                COUNT(*) AS sequence_count
            FROM meta_data md
            JOIN m49_country m on m.m49_code = md.country_validated
            {where_sql}
            GROUP BY
                md.country_validated,
                m.display_name
            ORDER BY sequence_count ASC;"""
        print(query, params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, params)
            # metadata_countries = cursor.fetchall()
            metadata_countries = dictfetchall(cursor)


        return metadata_countries

    def _add_filters(self):

        comparison_filters = {
                'length_lower': ('real_length', '>='),
                'length_upper': ('real_length', '<='),
                'collection_year_lower': ('collection_year', '>='),
                'collection_year_upper': ('collection_year', '<='),
                'creation_year_lower': ('create_date', '>='),
                'creation_year_upper': ('create_date', '<=')
            }

        where_clauses, params = [], []

        for key, value in self.filters.items():
            if key in comparison_filters:

                col, op = comparison_filters[key]
                where_clauses.append(f"{col} {op} %s")
                params.append(int(value))

            elif key == 'region':
                where_clauses.append(
                    "primary_accession IN ("
                    "SELECT m.primary_accession"
                    "FROM m49_country r "
                    "JOIN meta_data m ON m.country_validated = CAST(r.m49_code AS TEXT) "
                    "WHERE r.m49_region_id = %s)"
                )
                params.append(value)

            else:
                if isinstance(value, list):
                    placeholders = ', '.join(['%s'] * len(value))
                    where_clauses.append(f"{key} IN ({placeholders})")
                    params.extend(value)
                else:
                    where_clauses.append(f"{key} = %s")
                    params.append(value)

        where_str = ' AND '.join(where_clauses)

        print(where_str, params)

        return where_str, params

