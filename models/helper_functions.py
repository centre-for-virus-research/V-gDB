from collections import Counter

def dictfetchall(cursor):
    """Returns all rows from a cursor as a list of dictionaries."""
    columns = [col[0] for col in cursor.description]
    return [
        dict(zip(columns, row))
        for row in cursor.fetchall()
    ]


def build_clade_tree(clades):
    nodes = {}
    
    for clade in clades:
        major = clade['major_clade']
        minor = clade['minor_clade']
        if major not in nodes:
            nodes[major] = {
                'name':major,
                'text':major,
                'parent':None,
                'nodes':[]
            }
        if minor != None:
            nodes[major]['nodes'].append({
                'name':minor,
                'text':minor,
                'parent':major
            })
    
    tree = []
    for clade in nodes.values():
        if len(clade['nodes']) == 0:
            clade['nodes'] = None
    
    tree = list(nodes.values())

    return(tree)



def build_feature_tree(features):

    tree = []
    nodes = {}

    for feature in features:
        name = feature['name']
        if name == 'whole_genome,NULL':
            name = 'whole_genome'
        display_name = feature['description'] if feature['description'] is not None else name
        parent_name = feature['parent_name']
        nodes[name] = {
            'name': name,
            'text': display_name,
            'parent': parent_name,
            'nodes': []
        }
    for feature in features:
        name = feature['name']
        if name == 'whole_genome,NULL':
            name = 'whole_genome'
        parent_name = feature['parent_name']
        display_name = feature['description']

        if parent_name is None:
            tree.append(nodes[name])
        else:
            nodes[parent_name]['nodes'].append(nodes[name])

    def clean_tree(node):
        if len(node['nodes']) == 0:
            node['nodes'] = None
        else:
            for child in node['nodes']:
                clean_tree(child)

    for root in tree:
        clean_tree(root)
    

    return(tree)

def build_filter_query(filters):

    for key, value in filters:
        new_value = value.split(',') if ',' in value else [value]
        
        if key == 'gb_length':
            if new_value[0] and new_value[1]:
                where_clauses.append(f"gb_length BETWEEN %s AND %s")
                params.append(new_value[0], new_value[1])



                params.extend([new_value[0] or 0, new_value[1] or "(SELECT MAX(gb_length) FROM rabv_sequence)"])

        
        elif key in ['gb_create_date', 'gb_update_date']:
            if new_value[0] or new_value[1]:
                where_clauses.append(f"YEAR({key}) BETWEEN %s AND %s")
                params.extend([new_value[0] or 0, new_value[1] or f"(SELECT MAX(YEAR({key})) FROM rabv_sequence)"])

        elif key == 'collection_year':
            if (new_value[0] or new_value[1]):
                if new_value[0]:
                    where_clauses.append("earliest_collection_year >= %s")
                    params.append(new_value[0])
                if new_value[1]:
                    where_clauses.append("latest_collection_year <= %s")
                    params.append(new_value[1])
        elif key in ['m49_region_id', 'm49_sub_region_id', 'development_status'] and new_value[0]:
                    country_clauses.append(f"{key} IN %s")
                    country_params.append(new_value)

        elif key in ['is_ldc', 'is_lldc', 'is_sids'] and value:
            country_clauses.append(f"{key} = %s")
            country_params.append(True if value == 'true' else False)

        elif value and key in ['major_clade', 'minor_clade', 'm49_country_id']:
            where_clauses.append(f"{key} IN %s")
            params.append(new_value)

        elif value:
            where_clauses.append(f"{key} LIKE %s")
            params.append(value)


