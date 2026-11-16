from django.db import connections
from models.helpers import * 
from Bio.Seq import Seq
from Bio import Phylo
from .maddog.plot_reference_tree import build_maddog_tree

def _build_lineages_tree(clades):
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

class Lineages:
    """
    A class to analyze genetic mutations by comparing sequence alignments
    to a master reference sequence in a given database.
    """


    def __init__(self, database):
        """
        Initialize the Mutations object.

        Args:
            database (str): Django database alias to use for queries.
            reference_sequence (str, optional): Accession ID of the reference sequence. Defaults to 'NC_001542'.
        """
        self.database = database
    
    def get_lineages(self):

        query = "SELECT * FROM genotypes ORDER BY major_clade ASC, minor_clade ASC;"

        with connections[self.database].cursor() as cursor:
            cursor.execute(query)
            lineages = dictfetchall(cursor)

            lineages_tree = _build_lineages_tree(lineages)

        return lineages_tree

    def get_maddog_lineages(self):

        query = 'SELECT newick FROM trees WHERE name="full_tree_midpoint_rooted"'
      

        with connections[self.database].cursor() as cursor:
            cursor.execute(query)
            maddog_tree = dictfetchall(cursor)

            maddog_tree = build_maddog_tree(maddog_tree[0]["newick"], self.database)

        return maddog_tree


