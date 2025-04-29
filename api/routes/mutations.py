from rest_framework.decorators import api_view
from rest_framework.response import Response


from models.mutations import Mutations


@api_view(['GET'])
def get_mutations(request):

    database = request.headers.get('database', 'RABV_NEW')
    database = 'RABV_NEW'

    params = dict(request.GET.items())

    hosts = params["host"].split(',')
    region = params["region"]
    codons = params["codon"].split(',')

    mutations = Mutations(database=database,
                                hosts=hosts, 
                                codons=codons,
                                region=region 
                                )
    data = mutations.get_mutations()

    return Response(data)

@api_view(['GET'])
def get_mutation_regions_and_codons(request):

    database = request.headers.get('database', 'RABV_NEW')
    mutations = Mutations(database=database)
    data = mutations.get_mutation_regions_and_codons()
   
    return Response(data)

@api_view(['GET'])
def get_sequence_mutation(request, primary_accession):
    database = request.headers.get('database', 'RABV_NEW')

    # get alignment and ref sequence
