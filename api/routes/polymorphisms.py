from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status
from django.db import connections
from django.http import HttpResponse
import datetime
import os

from models.helpers import *
from models.polymorphisms import Polymorphisms

@api_view(['GET'])
def get_polymorphisms(request):

    database = request.headers.get('database', 'default')
    
    params = dict(request.GET.items())
    polymorphisms = Polymorphisms(database=database)
    try:
        data = polymorphisms.get_polymorphisms(params)
    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )
        
    return Response(data)


@api_view(['GET'])
def get_polymorphism(request, id):

    database = request.headers.get('database', 'default')
    

    polymorphisms = Polymorphisms(database=database)
    try:
        data = polymorphisms.get_polymorphism(id)
    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )
        
    return Response(data)