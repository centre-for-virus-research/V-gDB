from django.db import connections
from models.helpers import dictfetchall  # Ensure this function is imported

class Features:
    def __init__(self, database):
        self.database = database  

    def get_features(self):

        with connections[self.database].cursor() as cursor:
            # cursor.execute('SELECT name FROM genes WHERE parent_name IS NOT NULL;')
            cursor.execute('SELECT product FROM features where accession=%s', ['NC_001542'])
            features = dictfetchall(cursor)

        return features
    
    def get_feature(self, primary_accession):
        with connections[self.database].cursor() as cursor:
            cursor.execute('SELECT * FROM features WHERE accession = %s', [primary_accession])

            features = dictfetchall(cursor)

        return features


    def get_clades_tree(self):
        # with connections[self.database].cursor() as cursor:
        #     cursor.execute('SELECT * FROM genotypes')

        #     features = dictfetchall(cursor)

        # nodes = {}

        # for clade in features:
        #     major = clade['major_clade']
        #     minor = clade['minor_clade']
        #     if major != "NULL":
        #         if major not in nodes: 
        #             nodes[major] = {
        #                 'name':major,
        #                 'text':major,
        #                 'parent':None,
        #                 'nodes':[]
        #             }
        #         if minor != None:
        #             nodes[major]['nodes'].append({
        #                 'name':minor,
        #                 'text':minor,
        #                 'parent':major
        #             })
        
        tree = []
        # for clade in nodes.values():
        #     if len(clade['nodes']) == 0:
        #         clade['nodes'] = None
        
        # tree = list(nodes.values())
        return tree
    
