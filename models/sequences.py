from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

    
from models import sequences_helpers as sh
from models import polymorphisms_helpers as ph

from models.filter_helper import FilterHelper

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

    def __build_pagination_query(self, where_clauses, filter_params, next_cursor, prev_cursor, items_per_page):

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
        
        return query, pagination_params


    def get_sequences(self, next_cursor, prev_cursor, items_per_page):
        
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
        print(query, pagination_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, pagination_params)
            results = dictfetchall(cursor)

        if prev_cursor:
            results.reverse()

        total_count = self.__get_sequences_count(filter_where_sql, filter_params)

        next_cursor_new, prev_cursor_new = None, None
        if results: 
            next_cursor_new = results[-1]["primary_accession"]
            prev_cursor_new = results[0]["primary_accession"]

        results_dict = {
                            "data": results,
                            "total_count": total_count,
                            "next_cursor": next_cursor_new,
                            "prev_cursor": prev_cursor_new,
                        }
        
        return results_dict

    def get_sequences_meta_data_download(self):
        
        where_clauses, filter_params = [], []
        columns_str = '*'

        ############# BUILD WHERE string if there are filters
        if self.filters:
            if ("metadata_columns" in self.filters.keys()):
                columns_str = self.filters["metadata_columns"]
                del self.filters["metadata_columns"]
            if self.filters:
                # where_str, filter_params = self._add_filters_refactored()
                filterHelper = FilterHelper(filters = self.filters, database=self.database)
                where_str, filter_params = filterHelper.add_filters()
                where_clauses.append(where_str)
        filter_where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""
        #############

        query = f"""
                    SELECT {columns_str}
                    FROM meta_data
                    {filter_where_sql}
                """
        print(query, filter_params)

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, filter_params)
            results = dictfetchall(cursor)

        return results
    
    def get_sequences_download(self):
        
        where_clauses, filter_params = [], []

        ############# BUILD WHERE string if there are filters
        if self.filters:
            # where_str, filter_params = self._add_filters_refactored()
            filterHelper = FilterHelper(filters = self.filters, database=self.database)
            where_str, filter_params = filterHelper.add_filters()
            where_clauses.append(where_str)
        filter_where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""
        #############

        query = f"""
                    SELECT header as sequence_id, sequence as alignment
                    FROM sequences
                    WHERE header IN (
                                        SELECT primary_accession
                                        FROM meta_data
                                        {filter_where_sql}
                                    )
                """

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, filter_params)
            results = dictfetchall(cursor)

        sequences = self.__build_alignment_file(results)

        return sequences


    def get_sequences_alignment(self, start_coordinate, end_coordinate, sequence_type, product):

        where_clauses, filter_params = [], []

        if self.filters:
            # where_str, filter_params = self._add_filters_refactored()
            filterHelper = FilterHelper(filters = self.filters, database=self.database)
            where_str, filter_params = filterHelper.add_filters()
            where_clauses.append(where_str)

        filter_where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

        filtered_primary_accession_query = f"""
                                                SELECT primary_accession 
                                                FROM meta_data
                                                {filter_where_sql}
                                            """
        with connections[self.database].cursor() as cursor:

            if product: 
                alignment_query = f"""
                    SELECT sa.*, f.*
                    FROM sequence_alignment sa
                    JOIN features f ON f.accession = sa.sequence_id
                    WHERE sa.sequence_id IN ({filtered_primary_accession_query}) AND f.product = %s
                    """
                filter_params.append(product)
            else:
                alignment_query = f"""
                        SELECT sa.*
                        FROM sequence_alignment sa
                        WHERE sa.sequence_id IN ({filtered_primary_accession_query})
                        """

            cursor.execute(alignment_query, filter_params)
            results = dictfetchall(cursor)

            if product:
                alignments = self.__parse_alignments(results, start_coordinate, end_coordinate, sequence_type)
            else: 
                alignments = self.__build_alignment_file(results)

        return alignments


        
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

            # Add hose data
            if meta_data["host_taxa_id"]:
                host_taxa_id = meta_data["host_taxa_id"]
                taxanomic_info = sh._get_taxa_from_host_taxa_id(cursor, host_taxa_id)

                if taxanomic_info:
                    result["taxanomic_info"] = taxanomic_info

            query_alignment_dict = sh._get_query_alignment_from_primary_accession(cursor, primary_accession)
        
            if query_alignment_dict:

                reference_accession = query_alignment_dict["alignment_name"]
                query_alignment_sequence = query_alignment_dict["alignment"]
                result["insertions"] = [query_alignment_dict["insertion"]]

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

            if self.database == 'HCV':
                mutations = ph._get_polymorphim_sequence(cursor, primary_accession)
                if mutations:
                    result["mutations"] = mutations

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

        where_clauses = []
        params = []

        if self.filters:
            # where_str, filter_params = self._add_filters_refactored()
            filterHelper = FilterHelper(filters = self.filters, database=self.database)
            where_str, filter_params = filterHelper.add_filters()
            where_clauses.append(where_str)
            params.extend(filter_params)
        
        # Build query
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

                    
        query = f"""
                SELECT md.country_validated as m49_code,
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
            metadata_countries = dictfetchall(cursor)

        return metadata_countries


    # def _add_filters_refactored(self):
    #     comparison_filters = {
    #             'length_lower': ('real_length', '<='),
    #             'length_upper': ('real_length', '>='),
    #             'collection_year_lower': ('collection_year', '<='),
    #             'collection_year_upper': ('collection_year', '>='),
    #             'creation_year_lower': ('create_date', '>='),
    #             'creation_year_upper': ('create_date', '<=')
    #         }
        
    #     taxonomy_filters = ['phylum', 'class', 'order_category', 'family', 'genus', 'species', 'host']
    #     standard_filters = ['primary_accession', 'isolate', 'exclusion_status', 'accession_type', 'country_validated']
    #     clade_filters = ['EPA_major_clade', 'EPA_minor_clade']
    #     region_filters = ['country', 'm49_region_id', 'm49_sub_region_id', 'm49_intermediate_region_id', 'm49_code']
    #     genome_coverage_filters = ['full_genome', 'nucleoprotein', 'phosphoprotein', 'm2_protein', 'glycoprotein', 'l_protein', 'coreprotein', 'envelope_protein_E1', 'envelope_protein_E2', 'protein_p7', 'NS2', 'NS3', 'NS4A', 'NS4B', 'NS5B']


    #     where_clauses, params = [], []
    #     taxonomy_clauses, taxonomy_params = [], []
    #     region_clauses, region_params = [], []
    #     genome_coverage_clauses, genome_coverage_params = [], []

    #     for key, value in self.filters.items():
    #         exclude = False
    #         if key.startswith("exclude_"): 
    #             exclude = True
    #             key=key[8:]
    #             if key == 'taxa':
    #                 continue
    #             if key == 'clades':
    #                 continue
    #             if key == 'coverage':
    #                 continue

    #         if key in comparison_filters:
    #             col, op = comparison_filters[key]
    #             where_clauses.append(f"{col} {op} %s")
    #             params.append(int(value))

    #         elif key in taxonomy_filters:
    #             clause, param = sh._add_standard_filters(key, value, exclude)
    #             taxonomy_clauses.append(clause)
    #             taxonomy_params.extend(param)

    #         elif key in region_filters:
    #             if key=='country':
    #                 key='display_name'
    #             clause, param = sh._add_standard_filters(key, value, exclude)
    #             region_clauses.append(clause)
    #             region_params.extend(param)

    #         elif key in genome_coverage_filters:
    #             if key == 'full_genome':
    #                 clause, param = sh._add_standard_filters('calculated_genome_coverage', value, exclude)
    #                 where_clauses.append(clause)
    #                 params.extend(param)
    #             else:
    #                 protein_map = {
    #                     'nucleoprotein': 'nucleoprotein N',
    #                     'phosphoprotein': 'phosphoprotein M1',
    #                     'm2_protein': 'M2 protein',
    #                     'glycoprotein': 'transmembrane glycoprotein G',
    #                     'l_protein': 'L protein',
    #                     'coreprotein': 'core protein',
    #                     'envelope_protein_E1': 'envelope protein E1',
    #                     'envelope_protein_E2': 'envelope protein E2',
    #                     'protein_p7': 'protein p7',
    #                     'NS2': 'nonstructural protein NS2',
    #                     'NS3': 'protease/helicase protein NS3',
    #                     'NS4A': 'nonstructural protein NS4A',
    #                     'NS4B': 'nonstructural protein NS4B',
    #                     'NS5B': 'RNA-dependent RNA polymerase NS5B'
    #                 }

    #                 product = protein_map.get(key)

    #                 if product:
    #                     clause, param = sh._add_genome_coverage_filters(product, value)
    #                     genome_coverage_clauses = clause
    #                     genome_coverage_params.extend(param)
                            

    #         elif key in clade_filters:
    #             if self.filters.get("exclude_clades"):
    #                 exclude = True
    #                 clause, param = sh._add_exclude_clade_filters(key, self.filters.get("EPA_major_clade"), self.filters.get("EPA_minor_clade"), exclude)
    #             else:
    #                 clause, param = sh._add_standard_filters(key, value, exclude)
    #             where_clauses.append(clause)
    #             params.extend(param)

    #         elif key in standard_filters:
    #             clause, param = sh._add_standard_filters(key, value, exclude)
    #             where_clauses.append(clause)
    #             params.extend(param)

    #     if genome_coverage_clauses:
    #         comparison = 'IN'
    #         if ("exclude_genome" in self.filters.keys()):
    #             comparison = 'NOT IN'
            
    #         clause = sh._add_genome_filter(genome_coverage_clauses, comparison)
    #         params.extend(genome_coverage_params)
                
    #         where_clauses.append(clause)
            

    #     if taxonomy_clauses:
            
    #         comparison = 'IN'
    #         if ("exclude_taxa" in self.filters.keys()):
    #             comparison = 'NOT IN'

    #         if ("species" in self.filters.keys()):
    #             clause, param = sh._add_taxonomy_species_filters(self.filters["species"], comparison)
    #             params.extend(param)
    #         elif ("host" in self.filters.keys()):
    #             clause, param = sh._add_taxonomy_host_filters(self.filters["host"], comparison)
    #             params.extend(param)
    #         else:
    #             clause = sh._add_taxonomy_filters(taxonomy_clauses, comparison)
    #             params.extend(taxonomy_params)
                
    #         where_clauses.append(clause)

    #     if region_clauses:
    #         comparison = 'IN'
    #         if ("exclude_region" in self.filters.keys()):
    #             comparison = 'NOT IN'
            
    #         clause = sh._add_region_filters(region_clauses, comparison)
    #         params.extend(region_params)
                
    #         where_clauses.append(clause)

    
    #     where_str = ' AND '.join(where_clauses)

    #     print(where_str, params)

    #     return where_str, params
                
    def __build_alignment_file(self, alignments):
        results = []
        for a in alignments:
            results.append(f">" + a["sequence_id"] + "\n" + a["alignment"] + "\n")


        return results

    def __parse_alignments(self, alignments, start_coordinate, end_coordinate, sequence_type):
        results= []

        for a in alignments:
            if not start_coordinate and not end_coordinate:  # Use the full region coordinates

                ref_start = int(a["cds_start"])
                ref_end = int(a["cds_end"])

                if sequence_type == "codon":  # User wants a codon
                    codon_start, codon_end = get_codon_labeling(ref_start, ref_end)
            else:  # User chooses the coordinates themselves

                ref_start = int(start_coordinate)
                ref_end = int(end_coordinate)

                if sequence_type == "codon":  # User wants a codon
                    ref_start = int(a["cds_start"])
                    ref_end = int(a["cds_end"])
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
    


    def __get_sequences_count(self, where_sql, params):

        query = f"""
                    SELECT COUNT(primary_accession)
                    FROM meta_data
                    {where_sql}
                """

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, params)
            total_count = cursor.fetchone()[0]

        return total_count