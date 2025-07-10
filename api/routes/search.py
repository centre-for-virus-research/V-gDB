from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections

from models.helpers import dictfetchall  # Ensure this function is imported


@api_view(['GET'])
def search_primary_accession_ids(request, query):

    database = request.headers.get('database', 'default')
    print(query)
    with connections[database].cursor() as cursor:

        cursor.execute("SELECT DISTINCT(primary_accession) FROM meta_data WHERE primary_accession LIKE %s;", [f"%{query}%"])
        data = dictfetchall(cursor)

    return Response(data)


@api_view(['GET'])
def search_isolate_ids(request, query):

    database = request.headers.get('database', 'default')

    with connections[database].cursor() as cursor:

        cursor.execute("SELECT DISTINCT(isolate) FROM meta_data WHERE isolate LIKE %s;", [f"%{query}%"])
        data = dictfetchall(cursor)

    return Response(data)


@api_view(['GET'])
def search_pubmed_ids(request, query):

    database = request.headers.get('database', 'default')

    with connections[database].cursor() as cursor:

        cursor.execute("SELECT DISTINCT(pubmed_id) FROM meta_data WHERE pubmed_id LIKE %s;", [f"%{query}%"])
        data = dictfetchall(cursor)

    return Response(data)


@api_view(['GET'])
def search_hosts(request, query):

    database = request.headers.get('database', 'default')

    with connections[database].cursor() as cursor:

        cursor.execute("SELECT DISTINCT(host) FROM meta_data WHERE host LIKE %s;", [f"%{query}%"])
        data = dictfetchall(cursor)

    return Response(data)

@api_view(['GET'])
def search_country(request, query):

    database = request.headers.get('database', 'default')

    with connections[database].cursor() as cursor:

        cursor.execute("SELECT DISTINCT(m49_code), display_name FROM m49_country WHERE display_name LIKE %s;", [f"%{query}%"])
        data = dictfetchall(cursor)

    return Response(data)

