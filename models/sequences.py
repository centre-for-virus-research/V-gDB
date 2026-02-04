from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

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
        # Base query (with filters if any)
        where_clauses = []
        params = []
        if self.filters:
            where_str, filter_params = self.__add_filters()
            where_clauses.append(where_str)
            params.extend(filter_params)
            # query = f"SELECT * FROM meta_data WHERE {where_str}"

            
        # Pagination logic
        order_by = "ORDER BY primary_accession"
        limit = "LIMIT %s"
        params.append(items_per_page)

        if next_cursor:
            where_clauses.append("primary_accession > %s")
            params.insert(-1, next_cursor)

        elif prev_cursor:
            order_by = "ORDER BY primary_accession DESC"
            if prev_cursor != "0":
                where_clauses.append("primary_accession < %s")
                params.insert(-1, prev_cursor)

        # Build query
        where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

                    
        query = f"""
                    WITH filtered AS (SELECT * FROM meta_data {where_sql})
                    SELECT *, (SELECT COUNT(primary_accession) FROM filtered) AS total_count
                    FROM filtered {order_by} {limit};
                """

        print(query, params)

        with connections["FLUV"].cursor() as cursor:
            cursor.execute(query, params)
            results = dictfetchall(cursor)

        # Normalize order for previous pagination
        if prev_cursor:
            results.reverse()

        print(results)

        return {
            "data": results,
            "total_count": results[0]["total_count"] if results else None,
            "next_cursor": results[-1]["primary_accession"] if results else None,
            "prev_cursor": results[0]["primary_accession"] if results else None,
        }

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
            meta_data = self.__get_meta_data_from_primary_accession(cursor, primary_accession)
            
            if not meta_data:
                raise ValueError(f"Sequence with primary_accession {primary_accession} not found")
            
            result["meta_data"] = meta_data

            sequence = self.__get_sequence_from_primary_accession(cursor, primary_accession)
            result.update(sequence)

            # Get regional info if country exists
            if meta_data["country_validated"]:
                country_code = int(meta_data["country_validated"])
                regions = self.__get_region_from_country_code(cursor, country_code)
                result["regions"] = regions

            # Get insertions
            insertions = self.__get_insertions_from_primary_accession(cursor, primary_accession)
            if insertions:
                result["insertions"] = insertions

            query_alignment_dict = self.__get_query_alignment_from_primary_accession(cursor, primary_accession)
        
            if query_alignment_dict:

                reference_accession = query_alignment_dict["alignment_name"]
                query_alignment_sequence = query_alignment_dict["alignment"]

                # Get aligned reference sequence
                reference_alignment_dict = self.__get_query_alignment_from_primary_accession(cursor, reference_accession)
                reference_alignment_sequence = reference_alignment_dict["alignment"]
                
                # Get features

                features = self.__get_features_from_primary_accession(cursor, primary_accession)

            
                alignment = {
                                "reference_accession": reference_accession,
                                "query_alignment_sequence":query_alignment_sequence,
                                "reference_alignment_sequence":reference_alignment_sequence,
                                "features": features,
                            }    
            
                result["alignment"] = alignment

        return result


    def get_map_data(self, data):
        filtered = [d["country"] for d in data if d.get("country") not in (None, "") ]
        print(filtered)
        with connections[self.database].cursor() as cursor:
            # cursor.execute(query, params)
            # metadata_countries = dictfetchall(cursor)
            # Get reference m49_country data
            cursor.execute('SELECT display_name, id, m49_code FROM m49_country')
            m49_country_data = cursor.fetchall()
        # Clean, merge and process country metadata with m49 reference data
        parsed_data = self.__parse_and_combine_country_data(filtered, m49_country_data)
        sequence_counts = self.__count_sequence_occurance_by_country(parsed_data)

        return sequence_counts

    def __add_filters(self):

        comparison_filters = {
                'length_lower': ('length', '>='),
                'length_upper': ('length', '<='),
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
                params.append(value)

            elif key == 'region':
                where_clauses.append(
                    "primary_accession IN ("
                    "SELECT m.primary_accession "
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

        return where_str, params


    def get_strains(self):
        # if not self.filters:
        #     with connections[self.database].cursor() as cursor:
        #         cursor.execute(' SELECT isolate, country, host, segment, collection_date, primary_accession FROM meta_data ORDER BY isolate, segment;')
        #         rows = cursor.fetchall()
        # # Structure data by isolate
        # data_by_isolate = defaultdict(lambda: {"host": None, "segments": {}})

        # for isolate, country, host, segment, collection_date, accession in rows:
        #     if data_by_isolate[isolate]["host"] is None:
        #         data_by_isolate[isolate]["host"] = host
        #         data_by_isolate[isolate]["country"] = country
        #         data_by_isolate[isolate]["collection_date"] = collection_date
        #     data_by_isolate[isolate]["segments"][segment] = accession

        # # Convert to a nice list/dict if needed
        # result = []
        # for isolate, data in data_by_isolate.items():
        #     entry = {
        #         "isolate": isolate,
        #         "country": data["country"],
        #         "collection_date": data["collection_date"],
        #         "host": data["host"],
        #         "segments": data["segments"]  # dict like {1: 'ACC001', 2: 'ACC002', ...}
        #     }
        #     result.append(entry)
        with connections[self.database].cursor() as cursor:
                cursor.execute(' SELECT * FROM isolates')
                result = dictfetchall(cursor)

        return result

    def get_strain(self, isolate):
        if isolate:
            with connections[self.database].cursor() as cursor:
                cursor.execute('SELECT * FROM meta_data where isolate=%s ORDER BY segment', [isolate])
                result = dictfetchall(cursor)
                # print(result)
                for l in result:
                    print(l["primary_accession"])
                    print("here0")
                    cursor.execute("SELECT * FROM sequence_alignment WHERE sequence_id=%s", [l["primary_accession"]])
                    alignment = dictfetchall(cursor)
                    print("here4")
                    if alignment:
                        l["alignment"] = alignment[0]
                        print("here")
                        # Add insertions
                        cursor.execute("SELECT * FROM insertions WHERE accession = %s", [l["primary_accession"]])
                        l["alignment"]["insertions"] = dictfetchall(cursor)
                        print("here2")
                        # Add features
                        print(l["alignment"]["alignment_name"])
                        cursor.execute("SELECT * FROM features WHERE accession=%s and reference_accession = %s ORDER BY cds_start", [l["primary_accession"], l["alignment"]["alignment_name"]])
                        l["alignment"]["features"] = dictfetchall(cursor)
                        print("here3")

                        # Add reference sequence
                        cursor.execute("SELECT alignment FROM sequence_alignment WHERE sequence_id = %s", [l["alignment"]["alignment_name"]])
                        l["alignment"]["ref_seq"] = dictfetchall(cursor)[0]["alignment"]

                # print(result)
                
        
        
        return result

    
    def __get_aligned_sequences_and_features_from_reference(cursor, reference_accession):

        query = f"""
                    SELECT * 
                    FROM sequence_alignment sa
                    JOIN features f ON f.accession = sa.sequence_alignment
                    WHERE sa.alignment_name = %s
                """
        params = [reference_accession]
        results = fetch_all(cursor, query, params)

        return results


    def get_reference_sequence(self, primary_accession):

        result = {}

        # Grabbing reference features
        with connections[self.database].cursor() as cursor:

            features = self.__get_features_from_primary_accession(cursor, primary_accession)
            # cursor.execute("SELECT * FROM features WHERE accession = %s", [primary_accession])

            # features = dictfetchall(cursor)
            print(features)
        
        if not features:
            raise ValueError("Reference sequence with primary_accession {primary_accession} not found")
        
        
        for feature in features:
            codons = get_codon_labeling(feature["cds_start"], feature["cds_end"])
            feature["codon_start"] = codons[0]
            feature["codon_end"] = codons[1]
            
        result["features"] = features

        with connections[self.database].cursor() as cursor:

            aligned_sequences = self.__get_aligned_sequences_and_features_from_reference
            
            



            cursor.execute("select * from sequence_alignment where alignment_name = %s", [primary_accession])
            result["aligned_sequences"] = dictfetchall(cursor)

            for i in range(0,len(result["aligned_sequences"])):
                cursor.execute("SELECT * FROM features WHERE accession = %s", [result["aligned_sequences"][i]["sequence_id"]])
                result["aligned_sequences"][i]["features"] = dictfetchall(cursor)


            # cursor.execute("SELECT sequence FROM sequences WHERE header=%s", [primary_accession])
            cursor.execute("SELECT alignment from sequence_alignment where sequence_id=%s", [primary_accession])
            result["ref_sequence"] = dictfetchall(cursor)[0]["alignment"]
            result["aligned_sequences"]

        return result

    def get_reference_sequence_old(self, primary_accession):

        result = {}

        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT * FROM features WHERE accession = %s", [primary_accession])

            features = dictfetchall(cursor)
        
        if not features:
            raise ValueError("Reference sequence with primary_accession {primary_accession} not found")
        
        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT length FROM meta_data WHERE primary_accession = %s", [primary_accession])
            max_whole_genome = cursor.fetchone()[0]

        
        for feature in features:
            codons = get_codon_labeling(feature["cds_start"], feature["cds_end"])
            feature["codon_start"] = codons[0]
            feature["codon_end"] = codons[1]
            
        
        # features.append({'ref_seq_name':primary_accession, 'feature_name':None, "ref_start":1, "ref_end":max_whole_genome, "product":"Whole Genome", "protein_id":None})
        result["features"] = features

        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT * FROM genes")
            result["genes"] = dictfetchall(cursor)

            cursor.execute("SELECT * FROM sequence_alignment WHERE alignment_name = %s", [primary_accession])
            result["aligned_sequences"] = dictfetchall(cursor)

            # cursor.execute("SELECT sequence FROM sequences WHERE header=%s", [primary_accession])
            cursor.execute("SELECT alignment from sequence_alignment where sequence_id=%s", [primary_accession])
            result["sequence"] = dictfetchall(cursor)[0]["alignment"]

        return result



    def __parse_and_combine_country_data(self, meta_data, m49_data):
        """
        Matches country names from metadata with m49_country data and enriches entries.

        Parameters:
            meta_data: List of country records from meta_data table.
            m49_data: List of tuples (display_name, id, m49_code) from m49_country.

        Returns:
            Enriched meta_data entries with m49_code and country id.
        """
        # Build lookup table from m49_country
        m49_lookup = {
            key: (m49_code, country_id)
            for display_name, country_id, m49_code in m49_data
            for key in (display_name, country_id)
        }

        for entry in meta_data:
            country = entry
            if country:
                parsed_country = 'Vietnam' if country.split(":")[0] == 'Viet Nam' else country.split(":")[0]
                entry["parsed_country"] = parsed_country
                entry["m49_code"], entry["id"] = m49_lookup.get(parsed_country, (None, None))

        return meta_data

    def __count_sequence_occurance_by_country(self, meta_data):
        """
        Counts the number of sequences per m49 code and returns enriched, unique country entries.

        Parameters:
            meta_data: List of enriched meta_data entries.

        Returns:
            List of unique countries with sequence counts.
        """
        # Count how many times each m49_code appears
        print(meta_data)
        m49_codes = [entry["m49_code"] for entry in meta_data]
        m49_code_counts = Counter(m49_codes)

        seen_m49_codes = set()
        unique_data = []

        for entry in meta_data:
            entry["sequence_count"] = m49_code_counts[entry["m49_code"]]
            if entry["m49_code"] not in seen_m49_codes:
                unique_data.append(entry)
                seen_m49_codes.add(entry["m49_code"])

        return unique_data


    def __get_meta_data_from_primary_accession(self, cursor, primary_accession):

        query = "SELECT * FROM meta_data WHERE primary_accession = %s;"
        params = [primary_accession]
        meta_data = fetch_one(cursor, query, params)
        
        return meta_data
    
    def __get_sequence_from_primary_accession(self, cursor, primary_accession):

        query = "SELECT sequence FROM sequences WHERE header=%s"
        params = [primary_accession]

        sequence = fetch_one(cursor, query, params)
        return sequence
    
    def __get_features_from_primary_accession(self, cursor, primary_accession):
        query = f"""
                    SELECT * 
                    FROM features 
                    WHERE accession=%s 
                    ORDER BY cds_start
                """
        params = [primary_accession]
        features = fetch_all(cursor, query, params)

        return features
    
    def __get_insertions_from_primary_accession(self, cursor, primary_accession):
        query = "SELECT * FROM insertions WHERE accession = %s"
        params = [primary_accession]

        insertions = fetch_all(cursor, query, params)
        return insertions
    
    def __get_query_alignment_from_primary_accession(self, cursor, primary_accession):

        query = "SELECT * FROM sequence_alignment WHERE sequence_id=%s"
        params = [primary_accession]

        query_alignment = fetch_one(cursor, query, params)
        return query_alignment

    def __get_region_from_country_code(self, cursor, country_code):

        query = "SELECT * FROM m49_country WHERE m49_code=%s"
        params = [country_code]

        regions = fetch_one(cursor, query, params)
        return regions