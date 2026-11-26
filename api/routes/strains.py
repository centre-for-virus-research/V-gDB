from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status
from django.db import connections
from django.http import HttpResponse
import datetime
import os

from models.helpers import *
from models.strains import Strains
from django.http import StreamingHttpResponse

@api_view(['GET'])
def get_strains(request):

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

    strains = Strains(database=database, filters=params)
    try:
        data = strains.get_strains(next_cursor, prev_cursor, items_per_page)
    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )
        
    return Response(data)


@api_view(['GET'])
def get_strain(request):

    database = request.headers.get('database', 'default')

    params = dict(request.GET.items())

    parsed_strain = params["strain_id"]
    strains = Strains(database=database)



    try:
        data = strains.get_strain(parsed_strain)
    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )

    return Response(data)

