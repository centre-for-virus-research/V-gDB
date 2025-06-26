from django.db import connections
from collections import Counter
from models.helpers import dictfetchall 

class Statistics:
    """
    A class to retrieve and process statistical information from a specified database.
    """

    def __init__(self, database, filters=None):
        """
        Initialize with the target database name.
        """
        self.database = database
        self.filters = filters  

    def get_global_distribution_of_sequences(self):
        """
        Returns a list of unique countries (with m49 codes) and the number of sequences per country.
        """

        if not self.filters:
           with connections[self.database].cursor() as cursor:
                # Get country values from metadata
                cursor.execute('SELECT country FROM meta_data WHERE country IS NOT NULL')
                metadata_countries = dictfetchall(cursor)

                # Get reference m49_country data
                cursor.execute('SELECT display_name, id, m49_code FROM m49_country')
                m49_country_data = cursor.fetchall()
        else:
            where_clauses = []
            params = []

            def add_filter_clause(key, value, operator='='):
                """Add WHERE clause for numeric or comparison-based filters."""
                where_clauses.append(f"{key} {operator} %s")
                params.append(value)

            def add_filter_in_clause(key, value):
                """Add simple equality filter for strings/other fields."""
                where_clauses.append(f"{key} = %s")
                params.append(value)

            for key, value in self.filters.items():
                if key == 'length_lower':
                    add_filter_clause('length', value, operator='>=')
                elif key == 'length_upper':
                    add_filter_clause('length', value, operator='<=')
                else:
                    add_filter_in_clause(key, value)

            where_str = ' AND '.join(where_clauses)
            query = f"SELECT country FROM meta_data WHERE {where_str} AND country IS NOT NULL"
            print(query, params)
            with connections[self.database].cursor() as cursor:
                cursor.execute(query, params)
                metadata_countries = dictfetchall(cursor)
                # Get reference m49_country data
                cursor.execute('SELECT display_name, id, m49_code FROM m49_country')
                m49_country_data = cursor.fetchall()
                print(m49_country_data)

            

        # Clean, merge and process country metadata with m49 reference data
        parsed_data = self.__parse_and_combine_country_data(metadata_countries, m49_country_data)
        sequence_counts = self.__count_sequence_occurance_by_country(parsed_data)

        return sequence_counts

    def get_sequences_count(self):
        """
        Returns the total number of sequences in the database.
        """
        with connections[self.database].cursor() as cursor:
            cursor.execute('SELECT COUNT(primary_accession) as sequences_count FROM meta_data;')
            sequences_count = cursor.fetchone()

        return sequences_count[0]

    def get_reference_sequences_count(self):
        """
        Returns the number of sequences that are also reference alignments.
        """
        with connections[self.database].cursor() as cursor:
            cursor.execute("""
                SELECT count(primary_accession) 
                FROM meta_data 
                WHERE primary_accession IN (SELECT alignment_name FROM sequence)
            """)
            sequences_count = cursor.fetchone()

        return sequences_count[0]

    def get_max_min_sequence_length(self):
        """
        Returns the minimum and maximum sequence lengths from the database.
        """
        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT MIN(length), MAX(length) FROM meta_data;")
            max_min = cursor.fetchone()
            
        return max_min

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
            country = entry.get("country")
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
