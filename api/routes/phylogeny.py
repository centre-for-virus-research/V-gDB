from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections

from Bio import Phylo
from io import StringIO
from models.helpers import *
from urllib.parse import quote
@api_view(['GET'])
def get_tree(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    print(database)
    with open("/Volumes/My Passport/CVR/gdb/Flu/trees/seg1_cluster_rep.nwk", "r") as f:
        data = f.read()

    # print(data)
    # tree = Phylo.read(StringIO(data), "newick")

    # # Leaf labels
    # leaf_labels = [term.name for term in tree.get_terminals()]
    # # leaf_labels = [quote(name) for name in leaf_labels]
    # # print("Leaf labels:", leaf_labels)
    # placeholders = ",".join("?" for _ in leaf_labels)
    # placeholders = ', '.join(['%s'] * len(leaf_labels))
    # # leaf_labels_quoted = [f'"{x}"' for x in leaf_labels]
    # # placeholders = ",".join(leaf_labels)
    # sql = f"""
    #     SELECT Parsed_strain, host
    #     FROM meta_data
    #     WHERE Parsed_strain IN ({placeholders})
    # """

    # print(sql)
    # print(leaf_labels[0])
    

    # with connections[database].cursor() as cursor:
    #     cursor.execute(sql, leaf_labels)
    #     result = dictfetchall(cursor)
    
    # print(result)
    result = None
    # data = {"tree":data, "meta_data":result}
    data = {"tree":data}
        
    return Response(data)
