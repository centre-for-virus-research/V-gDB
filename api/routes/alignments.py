from rest_framework.decorators import api_view
from rest_framework.response import Response

from models.helpers import *
from django.http import HttpResponse
from models.helpers import *

from urllib.parse import unquote

# from models.alignment import Alignment
from models.sequences import Sequences


import datetime
import os
import json


@api_view(['GET'])
def download_alignments(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())

    region = None
    sequence_type = None
    start_coordinate = None
    end_coordinate = None

    if "region" in params:
        region = params["region"] 
        del params["region"]

    if "sequenceType" in params:
        sequence_type = params["sequenceType"]
        del params["sequenceType"]

    if "startCoordinate" in params:
        start_coordinate = params["startCoordinate"]
        del params["startCoordinate"]

    if "endCoordinate" in params:
        end_coordinate = params["endCoordinate"]
        del params["endCoordinate"]
    if "items_per_page" in params:
        del params["items_per_page"]
    if "next_cursor" in params:
        del params["next_cursor"]
    if "prev_cursor" in params:
        del params["prev_cursor"]
    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)
    data = sequences.get_sequences_alignment(start_coordinate, end_coordinate, sequence_type, region)

    

    file_name = str(datetime.datetime.now().strftime('%Y-%m-%d')) + '_alignments.fasta'
    build_fasta_file(data, file_name)

    with open(file_name, 'r') as file:
        response = HttpResponse(file, content_type='text')
        response['Content-Disposition'] = f'attachment; filename="{file_name}"'
        os.remove(file_name)

    return response


