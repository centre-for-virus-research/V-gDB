from rest_framework.decorators import api_view
from rest_framework.response import Response
from models.mutations import Mutations

@api_view(['GET'])
def get_mutations2(request):

    database = request.headers.get('database')
    params = dict(request.GET.items())

    hosts = params["host"].split(',')
    region = params["region"]
    codons = params["codon"].split(',')

    mutations = Mutations(database=database)
    data = mutations.get_mutations(hosts, codons, region)

    return Response(data)

@api_view(['GET'])
def get_mutations(request):

    database = request.headers.get('database')
    params = dict(request.GET.items())

    sequence_ids_param = params.get("sequence_ids")  # returns None if not present
    sequence_ids = sequence_ids_param.split('') if sequence_ids_param else None

    host_params = params.get("host") 
    hosts = host_params.split(',') if host_params else None

    include_metadata_param = params.get("include_metadata")  # returns None if not present
    include_metadata = include_metadata_param.split('') if include_metadata_param else True


    region = params["region"]
    codons = params["codon"].split(',')


    mutations = Mutations(database=database)
    data = mutations.get_mutations(codons=codons, region=region, include_metadata=include_metadata, sequence_ids=sequence_ids, hosts=hosts)

    return Response(data)

@api_view(['GET'])
def get_mutation_regions_and_codons(request):

    database = request.headers.get('database')
    mutations = Mutations(database=database)
    data = mutations.get_mutation_regions_and_codons()
   
    return Response(data)

