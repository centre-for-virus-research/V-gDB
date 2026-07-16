from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections

from models.phylogeny import Phylogeny



@api_view(['GET'])
def get_trees(request):

    database = request.headers.get('database', 'default')
    
    phylogeny = Phylogeny(database=database)

    tree = phylogeny.get_trees()
    
    return Response(tree)