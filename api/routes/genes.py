from rest_framework.decorators import api_view
from rest_framework.response import Response


from models.genes import Genes


@api_view(['GET'])
def get_genes_tree(request):

    database = request.headers.get('database', 'RABV_NEW')

    genes = Genes(database=database)
    data = genes.get_genes_tree()

    return Response(data)