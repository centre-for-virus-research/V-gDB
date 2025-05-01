from rest_framework.decorators import api_view
from rest_framework.response import Response

from models.helpers import *
from django.http import HttpResponse
from models.helpers import *

from urllib.parse import unquote

from models.alignment import Alignment
from models.sequences import Sequences


import datetime
import os
import json


@api_view(['GET'])
def download_alignments(request):

    database = request.headers.get('database', 'default')

    
    params = dict(request.GET.items())
    filters = json.loads(unquote(params["filters"])) if params["filters"] != 'undefined' else None
    sequences_helper = Sequences(database=database, filters=filters)
    data = sequences_helper.get_sequences_meta_data_by_filters()
    sequences = [d["primary_accession"] for d in data if "primary_accession" in d]
    region = params['region']
    nucleotide_or_codon = params['nucleotide_or_codon']
    start_coordinate = params['start']
    end_coordinate = params['end']


    alignment = Alignment(database=database,
                          sequences=sequences,
                        region=region, 
                        nucleotide_or_codon=nucleotide_or_codon, 
                        start_coordinate=start_coordinate, 
                        end_coordinate=end_coordinate)
    
    data = alignment.get_alignments()
    
    

    file_name = str(datetime.datetime.now().strftime('%Y-%m-%d')) + 'alignment.fasta'
    build_fasta_file(data, file_name)

    with open(file_name, 'r') as file:
        response = HttpResponse(file, content_type='text')
        response['Content-Disposition'] = 'attachment; filename='+file_name
        os.remove(file_name)

    return response


