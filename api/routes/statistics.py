from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections
from models.helper_functions import *
from django.http import HttpResponse
from models.helper_functions import *
from models.codon_labeling import get_kuiken2006_codon_labeling
# from Levenshtein import distance as levenshtein_distance
from difflib import SequenceMatcher

from models.statistics import Statistics

@api_view(['GET'])
def get_global_distribution_of_sequences(request):

    database = request.headers.get('database', 'default')

    statistics = Statistics(database=database)

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

    sequences_count = statistics.get_sequences_count()
    reference_sequences_count = statistics.get_reference_sequences_count()
    min_max_length = statistics.get_max_min_sequence_length

    data["sequences_count"] = sequences_count
    data["reference_sequences_count"] = reference_sequences_count
    data["min_length"] = min_max_length[0]
    data["max_length"] = min_max_length[1]

    return Response(data)


