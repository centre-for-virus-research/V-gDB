from django.db import connections
from collections import Counter
from models.helpers import dictfetchall  # Ensure this function is imported

class Statistics:
    def __init__(self, database):
        self.database = database  

    def get_global_distribution_of_sequences(self):

        with connections[self.database].cursor() as cursor:
            # Fetch all countries from meta_data
            cursor.execute('SELECT country FROM meta_data WHERE country IS NOT NULL')
            metadata_countries = dictfetchall(cursor)

            # Fetch all m49_country data in one query
            cursor.execute('SELECT display_name, id, m49_code FROM m49_country')
            m49_country_data = cursor.fetchall()

        parsed_data = self.__parse_and_combine_country_data(metadata_countries, m49_country_data)
        sequence_counts = self.__count_sequence_occurance_by_country(parsed_data)

        result = sequence_counts

        
        return result

    def get_sequences_count(self):

        with connections[self.database].cursor() as cursor:
            cursor.execute('SELECT COUNT(primary_accession) as sequences_count FROM meta_data;')
            sequences_count = cursor.fetchone()

        return sequences_count[0]

    def get_reference_sequences_count(self):
        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT count(primary_accession) FROM meta_data WHERE primary_accession IN (SELECT alignment_name FROM sequence)")
            sequences_count = cursor.fetchone()

        return sequences_count[0]

    def get_max_min_sequence_length(self):
        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT MIN(length), MAX(length) FROM meta_data;")
            max_min = cursor.fetchone()
            
        return max_min

    def __parse_and_combine_country_data(self, meta_data, m49_data):

        m49_lookup = {
            key: (m49_code, country_id)
            for display_name, country_id, m49_code in m49_data
            for key in (display_name, country_id)
        }

        # Process metadata countries
        for entry in meta_data:
            country = entry.get("country")
            if country:
                # Normalize and parse country name to match m49 format
                parsed_country = 'Vietnam' if country.split(":")[0] == 'Viet Nam' else country.split(":")[0]
                entry["parsed_country"] = parsed_country

                # Look up m49_code and id
                entry["m49_code"], entry["id"] = m49_lookup.get(parsed_country, (None, None))

        return meta_data



    def __count_sequence_occurance_by_country(self, meta_data):
        # Count occurrences of each m49_code
        m49_codes = [entry["m49_code"] for entry in meta_data]
        m49_code_counts = Counter(m49_codes)

        # Add sequence_count and remove duplicates
        seen_m49_codes = set()
        unique_data = []
        for entry in meta_data:
            entry["sequence_count"] = m49_code_counts[entry["m49_code"]]
            if entry["m49_code"] not in seen_m49_codes:
                unique_data.append(entry)
                seen_m49_codes.add(entry["m49_code"])

        return unique_data