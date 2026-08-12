from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections
from models.helpers import *
from django.http import HttpResponse
from models.helpers import *


from models.statistics import Statistics

@api_view(['GET'])
def get_global_distribution_of_sequences(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())

    for key, value in params.items():
        params[key] = value.split(',') if ',' in value else value

    statistics = Statistics(database=database, filters=params)

    try:
        data = statistics.get_global_distribution_of_sequences()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
        
    return Response(data)

@api_view(['GET'])
def get_statistics(request):

    database = request.headers.get('database', 'default')

    data = {}
    statistics = Statistics(database=database)
    pipeline_last_run = statistics.get_pipeline_last_run()
    sequences_count = statistics.get_sequences_count()
    # reference_sequences_count = statistics.get_reference_sequences_count()
    # recent_collection_year = statistics.get_recent_collection_year()
    # min_max_length = statistics.get_max_min_sequence_length
    data["pipeline_last_run"] = pipeline_last_run
    data["sequences_count"] = sequences_count
    # data["reference_sequences_count"] = reference_sequences_count
    # data["recent_collection_year"] = recent_collection_year
    # data["min_length"] = min_max_length[0]
    # data["max_length"] = min_max_length[1]

    return Response(data)