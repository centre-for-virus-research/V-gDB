from django.db import connections
import csv
from models.helper_functions import dictfetchall  # Ensure this function is imported

class Sequences:
    def __init__(self, database, primary_accession=None, filters=None):
        self.database = database  
        self.primary_accession = primary_accession
        self.filters = filters

    def get_sequences_meta_data(self):
        
        with connections[self.database].cursor() as cursor:
            cursor.execute('SELECT * FROM meta_data ORDER BY create_date DESC;')
            result = dictfetchall(cursor) 
        return result
    
    def get_sequence_meta_data(self):

        result = {}

        if not self.primary_accession:
            raise ValueError("Primary accession can not be blank")

        with connections[self.database].cursor() as cursor:

            cursor.execute("SELECT * FROM meta_data WHERE primary_accession = %s;", [self.primary_accession] )
            sequence = dictfetchall(cursor)
            
            if not sequence:
                raise ValueError("Sequence with primary_accession {} not found".format(self.primary_accession))
            
            else:
                result["meta_data"] = sequence[0]
                cursor.execute("SELECT * FROM sequence_alignment WHERE sequence_id=%s", [self.primary_accession])
                alignment = dictfetchall(cursor)
                
                if (alignment):
                    result["alignment"] = alignment[0]

                    cursor.execute("SELECT * FROM insertions WHERE accession = %s", [self.primary_accession])
                    result["alignment"]["insertions"] = dictfetchall(cursor)
                    
                    
                    cursor.execute("SELECT * \
                                FROM features \
                                WHERE accession = %s \
                                ORDER BY cds_start", ['EF437215'])
                    result["alignment"]["features"] = dictfetchall(cursor)

                    # cursor.execute("SELECT alignment FROM sequence_alignment WHERE sequence_id = %s", [result["alignment"]["alignment_name"]])
                    # result["alignment"]["ref_seq"] = dictfetchall(cursor)[0]["alignment"]

                    cursor.execute("SELECT sequence FROM sequences WHERE header = %s", [result["alignment"]["alignment_name"]])
                    result["alignment"]["ref_seq"] = dictfetchall(cursor)[0]["sequence"]

                if (result["meta_data"]["country"]):
                    cursor.execute("SELECT * FROM m49_country WHERE full_name=%s", [result["meta_data"]["country"].split(':')[0]])

                    country = dictfetchall(cursor)
                    if(country):
                        result["meta_data"]["region"] = country[0]

        return result

    def get_sequences_meta_data_by_filters(self):

        # If there are no filters inputted
        if not self.filters:
            result = self.get_sequences_meta_data()

        else:

            where_clauses = []
            params = []
            country_clauses = []
            country_params = []

            def add_filter_clause(key, value, operator='='):
                """Helper function to add a clause and its parameters to the query."""
                where_clauses.append(f"{key} {operator} %s")
                params.append(value)

            def add_filter_in_clause(key, value):
                where_clauses.append(f"{key} = %s")
                params.append(value)
                # where_clauses.append(f"{key} IN ({','.join(['%s'] * len(value))})")
                # params.extend(value)
            print(self.filters.items())
            for key, value in self.filters.items():

                if key == 'length_lower':
                    add_filter_clause('length', value, operator='>=')
                elif key == 'length_upper':
                    add_filter_clause('length', value, operator='<=')

                else:
                    add_filter_in_clause(key, value)

            where_str = ' AND '.join(where_clauses)
            query = f"SELECT * FROM meta_data WHERE {where_str} ORDER BY create_date DESC;"
            print(query, params)

            with connections[self.database].cursor() as cursor:
                cursor.execute(query, params)
                result = dictfetchall(cursor)

        return result


    def download_sequences_meta_data(self):

        with connections[self.database].cursor() as cursor:
            cursor.execute("SELECT * FROM meta_data")
            sequences = dictfetchall(cursor)

        return sequences


    def build_meta_data_csv_file(self, data):

        with open("tmp_file.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(data[0].keys())
            for i in data:
                w.writerow(i.values())
            f.close()


























    # database = request.headers.get('database', 'default')
    # where_clauses = []
    # params = []
    # country_clauses = []
    # country_params = []

    # print(dict(request.GET.items()))
    # params = dict(request.GET.items())
    # for key, value in params.items():
    #     print("HERE")
    #     params[key] = value.split(',') if ',' in value else value
    # print("DONE")
    # print(params)
    # for key, value in params.items():
    #     new_value = value
    #     # new_value = value.split(',') if ',' in value else value
    #     # new_value = value

    #     if key == 'gb_length_lower':
    #         where_clauses.append(f"length >= %s ")
    #         params.append(new_value)
    #     elif key == 'gb_length_upper':
    #         where_clauses.append(f"length <= %s ")
    #         params.append(new_value)

    #     elif key == 'gb_update_date_lower':
    #         where_clauses.append(f"YEAR(update_date) >= %s ")
    #         params.append(new_value)
    #     elif key == 'gb_update_date_upper':
    #         where_clauses.append(f"YEAR(update_date) <= %s ")
    #         params.append(new_value)

    #     elif key == 'creation_year_lower':
    #         where_clauses.append(f"YEAR(create_date) >= %s ")
    #         params.append(new_value)
    #     elif key == 'creation_year_upper':
    #         where_clauses.append(f"YEAR(create_date) <= %s ")
    #         params.append(new_value)

    #     elif key == 'collection_year_lower':
    #         where_clauses.append(f'STRFTIME("%Y", collection_date) >= %s ')
    #         params.append(new_value)
    #     elif key == 'collection_year_upper':
    #         where_clauses.append(f"YEAR(collection_date) <= %s ")
    #         params.append(new_value)

    #     elif key in ['primary_accession', 'host', 'isolate', 'pubmed_id']:
    #         new_value = [new_value] if isinstance(new_value, str) else new_value
    #         where_clauses.append(f"{key} IN ({','.join(["%s"] * len(new_value))})")
    #         params.extend(new_value)

    #     elif key in ['major_clade', 'minor_clade']:
    #         where_clauses.append(f"{key} IN %s")
            


    #     elif key in ['m49_region_id', 'm49_sub_region_id', 'm49_code', 'development_status']:
    #         country_clauses.append(f"{key} IN %s")
    #         country_params.append(new_value)

    #     elif key in ['is_ldc', 'is_lldc', 'is_sids']:
    #         country_clauses.append(f"{key} = %s")
    #         country_params.append(True if value == 'true' else False)


    # # Country-specific query
    # if country_clauses:
        
    #     country_query = f"SELECT DISTINCT(id) FROM m49_country WHERE {' AND '.join(country_clauses)}"
    #     print(country_query) 
    #     print(country_params)
    #     with connections[database].cursor() as cursor:
    #         cursor.execute(country_query, country_params)
    #         country_ids = [row[0] for row in cursor.fetchall()]
    #     if country_ids:
    #         where_clauses.append("m49_country_id IN %s")
    #         params.append(country_ids)

    # # Main query
    # where_str = ' AND '.join(where_clauses)
    # query = f"SELECT * FROM meta_data \
    #             WHERE {where_str} ORDER BY create_date DESC;" if where_clauses else "SELECT * FROM meta_data"
    # print(query) 
    # print(params)
    # with connections[database].cursor() as cursor:
    #     cursor.execute(query, params)
    #     result = dictfetchall(cursor)

        
    # return Response(result)