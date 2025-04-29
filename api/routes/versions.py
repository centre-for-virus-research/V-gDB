from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections
from models.helpers import *

from models.helpers import *


@api_view(['GET'])
def get_vgt_version(request):

    ########
    # List of meta data for all sequences 
    ########

    database = request.headers.get('database', 'default')

    with connections[database].cursor() as cursor:
        cursor.execute('SELECT * FROM project_settings where name LIKE "PROJECT_VERSION" or name LIKE "%EXTENSION_BUILD_DATE%";')
        result = dictfetchall(cursor)
        
    return Response(result)

@api_view(['GET'])
def get_meta_data_columns(request):
    database = request.headers.get('database', 'default')

    with connections[database].cursor() as cursor:
        cursor.execute('PRAGMA table_info(meta_data);')
        result = dictfetchall(cursor)
        
    return Response(result)

@api_view(['GET'])
def check_db_connection(request):
    try:
        # Try to connect to the database
        connections['default'].cursor()
        return Response(True)
    except:
        # If there's an error, the database is not connected
        return Response(False)