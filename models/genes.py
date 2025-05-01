from django.db import connections
import csv
from models.helpers import dictfetchall  # Ensure this function is imported

class Genes:
    def __init__(self, database):
        self.database = database  


    def get_genes_tree(self):
        with connections[self.database].cursor() as cursor:
            # cursor.execute("SELECT name, display_name, parent_name, description FROM genes;")
            cursor.execute('SELECT product FROM features where accession=%s', ['NC_001542'])
            result = dictfetchall(cursor)

        tree = self.__build_genes_tree(result)


        return tree

    def __build_genes_tree(self, genes):
        tree = []
        nodes = {}
        for gene in genes:
            name = gene['product']
            display_name = name.capitalize()
            # parent_name = gene['parent_name']
            nodes[name] = {
                'name': name,
                'text': display_name,
                'parent': None,
                'nodes': []
            }
        for gene in genes:
            name = gene['product']
            parent_name = None
            display_name = name

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

        return tree

    def __build_genes_tree_old(self, genes):

        tree = []
        nodes = {}

        for gene in genes:
            name = gene['name']
            if name == 'whole_genome,NULL':
                name = 'whole_genome'
            display_name = gene['description'] if gene['description'] is not None else name
            parent_name = gene['parent_name']
            nodes[name] = {
                'name': name,
                'text': display_name,
                'parent': parent_name,
                'nodes': []
            }
        for gene in genes:
            name = gene['name']
            if name == 'whole_genome,NULL':
                name = 'whole_genome'
            parent_name = gene['parent_name']
            display_name = gene['description']

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
        

        return tree