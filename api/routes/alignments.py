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

    database = request.headers.get('database', 'RABV_NEW')

    
    params = dict(request.GET.items())
    filters = json.loads(unquote(params["filters"]))

    if filters:
        sequences_helper = Sequences(database=database, filters=filters)
        data = sequences_helper.get_sequences_meta_data_by_filters()
        sequences = [d["primary_accession"] for d in data if "primary_accession" in d]

    else:
        sequences = params['sequences'].split(',')

    region = params['region']
    nucleotide_or_codon = params['nucleotide_or_codon']
    start_coordinate = params['start']
    end_coordinate = params['end']


    alignment = Alignment(database=database)
    
    data = alignment.get_alignments(sequences=sequences,
                                region=region, 
                                nucleotide_or_codon=nucleotide_or_codon, 
                                start_coordinate=start_coordinate, 
                                end_coordinate=end_coordinate)
    
    

    file_name = str(datetime.datetime.now().strftime('%Y-%m-%d')) + 'alignment.fasta'
    build_fasta_file(data, file_name)

    with open('my_fasta.fasta', 'r') as file:
        response = HttpResponse(file, content_type='text')
        response['Content-Disposition'] = 'attachment; filename='+file_name
        os.remove(file_name)

    return response




    
    

    # query = "/Users/dana/CVR/backend/models/mt862689.fa"
    # ref_seq = "/Users/dana/CVR/backend/models/ef437215.fa"
    # tmp_folder = "/Users/dana/CVR/backend/models/tmp/"
    # output_file = "query_tophits.tsv"
    # master_acc = "EF437215"
    # processor = BlastAlignment(query, ref_seq, tmp_folder, output_file, master_acc)
    # processor.run_makeblastdb()
    # processor.run_blastn()
    # # processor.process_non_segmented_virus()
    # print("done")

    # with open('/Users/dana/CVR/backend/models/tmp/query_tophits.tsv', 'r') as file:
    #     response = HttpResponse(file, content_type='text')
    #     response['Content-Disposition'] = 'attachment; filename=blast_results.tsv'
    #     # os.remove('my_fasta.fasta')
    return Response(job.id)

