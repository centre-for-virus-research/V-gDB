from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections
from django.http import HttpResponse
import datetime
import os

from models.helpers import *
from models.sequences import Sequences

@api_view(['GET'])
def get_sequences(request):

    database = request.headers.get('database', 'default')
    # sequences = Sequences(database=database)

    params = dict(request.GET.items())

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)
    try:
        data = sequences.get_sequences()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
        
    return Response(data)

@api_view(['GET'])
def get_sequences_meta_data(request):

    database = request.headers.get('database', 'default')
    # sequences = Sequences(database=database)

    params = dict(request.GET.items())

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)
    try:
        data = sequences.get_sequences_meta_data()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
        
    return Response(data)

@api_view(['GET'])
def get_sequences_alignment(request):

    database = request.headers.get('database', 'default')

    params = dict(request.GET.items())

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)
    try:
        data = sequences.get_sequences_alignment()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
        
    return Response(data)

@api_view(['GET'])
def get_reference_sequences_meta_data(request):

    database = request.headers.get('database', 'default')

    params = dict(request.GET.items())

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)
    try:
        data = sequences.get_reference_sequences_meta_data()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
        
    return Response(data)


@api_view(['GET'])
def get_sequence_meta_data(request, primary_accession):

    database = request.headers.get('database', 'default')
    sequences = Sequences(database=database)

    try:
        data = sequences.get_sequence_meta_data(primary_accession)
    except ValueError as e:
        print("HI")
        print(str(e))
        return Response({'message': str(e)}, status=404)

    return Response(data)


@api_view(['GET'])
def get_sequence_alignment(request, primary_accession):

    database = request.headers.get('database', 'default')
    sequences = Sequences(database=database)

    try:
        data = sequences.get_sequence_alignment(primary_accession)
    except ValueError as e:
        return Response({'message': str(e)}, status=404)

    return Response(data)


@api_view(['GET'])
def get_reference_sequence(request, primary_accession):

    if not primary_accession:
        return HttpResponse("Reference sequence not defined", status=404)

    database = request.headers.get('database', 'default')

    sequences = Sequences(database=database)
    data = sequences.get_reference_sequence(primary_accession)

    return Response(data)


@api_view(['GET'])
def get_map_metadata(request):

    database = request.headers.get('database', 'default')

    params = dict(request.GET.items())

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    sequences = Sequences(database=database, filters=params)
    try:
        data = sequences.get_sequences_meta_data()
        data = sequences.get_map_data(data)
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
        
    return Response(data)




@api_view(['GET'])
def download_sequences_meta_data(request):
    print("starting")
    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())

    print(params)

    sequences = Sequences(database=database, filters=params)
    if params:
        data = sequences.get_sequences_meta_data_by_filters()
    else:
        data = sequences.get_sequences_meta_data()
    
    
    print("ending")
    file_name = str(datetime.datetime.now().strftime('%Y-%m-%d')) + '_meta_data.csv'
    build_csv_file(data, file_name)
    print("THIS IS THE FOWNLOAD")
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


