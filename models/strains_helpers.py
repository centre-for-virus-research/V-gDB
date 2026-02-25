from .helpers import *

def _get_all_strains(cursor):

    query = "SELECT * FROM isolates;"
    data = fetch_all(cursor, query, params=None)
    return data


def _get_segment_meta_data_from_strain_id(cursor, strain_id):
    query = 'SELECT * FROM meta_data where Parsed_strain=%s ORDER BY segment'
    params = [strain_id]
    data = fetch_all(cursor, query, params)
    print(data)
    return data