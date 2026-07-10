from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status
from django.db import connections
from django.http import HttpResponse
import datetime
import os

from models.helpers import *
from models.sequences import Sequences

@api_view(['GET'])
def get_sequences(request):

    database = request.headers.get('database', 'default')
    prev_cursor = None
    next_cursor = None

    params = dict(request.GET.items())

    if "next_cursor" in params:
        next_cursor = params["next_cursor"]
        del params["next_cursor"]

    if "prev_cursor" in params:
        prev_cursor = params["prev_cursor"]
        del params["prev_cursor"]

    if "items_per_page" in params:
        items_per_page = params["items_per_page"]
        del params["items_per_page"]

    if "EPA_minor_clade" in params:
        if params["EPA_minor_clade"] == "null":
            del params["EPA_minor_clade"]

    if "country_validated" in params:
        if params["country_validated"][0] == "0":
            params["country_validated"] = params["country_validated"][1:]

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)
    try:
        data = sequences.get_sequences(next_cursor, prev_cursor, items_per_page)
    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )
        
    return Response(data)


@api_view(['GET'])
def get_sequence(request, primary_accession):

    database = request.headers.get('database', 'default')
    sequences = Sequences(database=database)

    try:
        data = sequences.get_sequence(primary_accession)
    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )

    return Response(data)


@api_view(['GET'])
def get_global_distribution_of_sequences(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())

    if "EPA_minor_clade" in params:
        if params["EPA_minor_clade"] == "null":
            del params["EPA_minor_clade"]
    
    if params:
        for key, value in params.items():
            print("PARAMS", key, value)
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)

    try:
        data = sequences.get_global_distribution_of_sequences()
    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )
        
    return Response(data)



@api_view(['GET'])
def get_reference_sequence(request, primary_accession):

    if not primary_accession:
        return HttpResponse("Reference sequence not defined", status=404)

    database = request.headers.get('database', 'default')

    sequences = Sequences(database=database)
    try:
        data = sequences.get_reference_sequence(primary_accession)
    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )
        

    return Response(data)



@api_view(['GET'])
def download_sequences_meta_data(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())

    if params:
        if "items_per_page" in params:
            del params["items_per_page"]

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)

    data = sequences.get_sequences_meta_data_download()
    
    
    file_name = str(datetime.datetime.now().strftime('%Y-%m-%d')) + '_meta_data.csv'
    build_csv_file(data, file_name)

    with open(file_name, 'r') as file:
        response = HttpResponse(file, content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename='+file_name
        os.remove(file_name)

        return response
    
@api_view(['GET'])
def download_sequences(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())

    if params:
        if "items_per_page" in params:
            del params["items_per_page"]

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)

    data = sequences.get_sequences_download()
    
    
    file_name = str(datetime.datetime.now().strftime('%Y-%m-%d')) + '_sequences.fasta'
    build_fasta_file(data, file_name)

    with open(file_name, 'r') as file:
        response = HttpResponse(file, content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename='+file_name
        os.remove(file_name)

        return response

    
# THIS IS BEING USED IN THE MUTATION GUI
# TODO: REMOVE this and use callback version to get host
@api_view(['GET'])
def get_host_species(request):
    database = request.headers.get('database', 'default')
    
    with connections[database].cursor() as cursor:
        cursor.execute('SELECT DISTINCT(host) FROM meta_data WHERE host IS NOT NULL;')
        result = dictfetchall(cursor)
    return Response(result)




@api_view(['GET'])
def advanced_filter(request, query):
    database = request.headers.get('database', 'default')
    query="SELECT * FROM meta_data WHERE"+query

    with connections[database].cursor() as cursor:
        cursor.execute(query)
        result = dictfetchall(cursor)
    return Response(result)


