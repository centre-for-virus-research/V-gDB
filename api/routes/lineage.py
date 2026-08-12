from rest_framework.decorators import api_view
from rest_framework.response import Response
from models.lineages import Lineages


@api_view(['GET'])
def get_lineage(request):

    database = request.headers.get('database', 'default')

    lineages = Lineages(database=database)
    data = lineages.get_lineages()

    return Response(data)