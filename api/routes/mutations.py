from rest_framework.decorators import api_view
from rest_framework.response import Response
from models.mutations import Mutations


@api_view(['GET'])
def get_adaptive_mutations(request):

    database = request.headers.get('database', 'default')

    mutations = Mutations(database=database)
    data = mutations.get_adaptive_mutations()

    return Response(data)

@api_view(['GET'])
def get_adaptive_mutations_chart(request, segment):

    database = request.headers.get('database', 'default')

    mutations = Mutations(database=database)
    data = mutations.get_adaptive_mutations_chart(segment)

    return Response(data)


@api_view(['GET'])
def get_adaptive_mutations_chart(request, segment):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    final_params = {}
    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value
        if "species" in params:
            final_params = {"species":params["species"]}
        if "genus" in params:
            final_params = {"genus":params["genus"]}
    print("params", final_params)

    mutations = Mutations(database=database)

    if database == "RABV" or "default" and params:
        data = mutations.get_adaptive_mutations_chart_RABV(final_params)
    else :
        data = mutations.get_adaptive_mutations_chart(segment)

    return Response(data)

@api_view(['GET'])
def get_mutations(request):

    database = request.headers.get('database', 'default')
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

    database = request.headers.get('database', 'default')
    mutations = Mutations(database=database)
    data = mutations.get_mutation_regions_and_codons()
   
    return Response(data)

