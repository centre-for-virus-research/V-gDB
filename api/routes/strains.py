from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections
from django.http import HttpResponse
import datetime
import os

from models.helpers import *
from models.strains import Strains

@api_view(['GET'])
def get_strains(request):

    database = request.headers.get('database', 'default')

    params = dict(request.GET.items())

    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value

    strains = Strains(database=database, filters=params)
    try:
        data = strains.get_strains()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
        
    return Response(data)

@api_view(['GET'])
def get_strain(request, strain_id):

    database = request.headers.get('database', 'default')
    strains = Strains(database=database)

    try:
        data = strains.get_strain(strain_id)
    except ValueError as e:
        print(str(e))
        return Response({'message': str(e)}, status=404)

    return Response(data)



