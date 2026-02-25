from rest_framework.decorators import api_view
from rest_framework.response import Response
from models.taxonomy import Taxonomy


@api_view(['GET'])
def get_phylum(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    if params:
        if "phylum" in params:
            del params["phylum"]
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    taxonomy = Taxonomy(database=database)
    data = taxonomy.get_phylum(params)

    return Response(data)

@api_view(['GET'])
def get_class(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    if params:
        if "class" in params:
            del params["class"]
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    taxonomy = Taxonomy(database=database)
    data = taxonomy.get_class(params)

    return Response(data)

@api_view(['GET'])
def get_order(request):
    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    if params:
        if "order_category" in params:
            del params["order_category"]
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    taxonomy = Taxonomy(database=database)
    data = taxonomy.get_order(params)

    return Response(data)

@api_view(['GET'])
def get_family(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    if params:
        if "family" in params:
            del params["family"]
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value


    taxonomy = Taxonomy(database=database)
    data = taxonomy.get_family(params)

    return Response(data)

@api_view(['GET'])
def get_genus(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    if params:
        if "genus" in params:
            del params["genus"]
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value


    taxonomy = Taxonomy(database=database)
    data = taxonomy.get_genus(params)

    return Response(data)

@api_view(['GET'])
def get_species(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    if params:
        if "species" in params:
            del params["species"]
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    taxonomy = Taxonomy(database=database)
    data = taxonomy.get_species(params)

    return Response(data)