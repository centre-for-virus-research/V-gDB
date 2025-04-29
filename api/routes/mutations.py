from rest_framework.decorators import api_view
from rest_framework.response import Response


from models.mutations import Mutations


@api_view(['GET'])
def get_mutations(request):

    database = request.headers.get('database')
    params = dict(request.GET.items())

    hosts = params["host"].split(',')
    region = params["region"]
    codons = params["codon"].split(',')

    mutations = Mutations(database=database)
    data = mutations.get_mutations(hosts, codons, region)

    return Response(data)

@api_view(['GET'])
def get_mutation_regions_and_codons(request):

    database = request.headers.get('database')
    mutations = Mutations(database=database)
    data = mutations.get_mutation_regions_and_codons()
   
    return Response(data)

