from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.db import connections
from models.helper_functions import *
from django.http import HttpResponse
from models.helper_functions import *
from models.codon_labeling import get_kuiken2006_codon_labeling
import os
import csv
# from Levenshtein import distance as levenshtein_distance
from difflib import SequenceMatcher
import json
from urllib.parse import unquote

from models.alignment import Alignment
from models.sequences import Sequences
# from models.alignment.blast_alignment import BlastAlignment

from models.tasks import Tasks

import django_rq

import rq

@api_view(['GET'])
def get_reference_sequences_meta_data(request):

    database = request.headers.get('database', 'default')
    alignment = Alignment(database=database)
    data = alignment.get_reference_sequences_meta_data()

    return Response(data)


@api_view(['GET'])
def get_reference_sequence(request, primary_accession):

    if not primary_accession:
        return HttpResponse("Reference sequence not defined", status=404)

    database = request.headers.get('database', 'default')

    alignment = Alignment(database=database)
    data = alignment.get_reference_sequence(primary_accession)

    return Response(data)


@api_view(['GET'])
def download_alignments(request):

    database = request.headers.get('database', 'RABV_NEW')
    database = 'RABV_NEW'

    
    params = dict(request.GET.items())
    filters = json.loads(unquote(params["filters"]))

    if filters:
        sequences_helper = Sequences(database=database, filters=filters)
        data = sequences_helper.get_sequences_meta_data_by_filters()
        sequences = [d["primary_accession"] for d in data if "primary_accession" in d]

    else:
        sequences = params['sequences'].split(',')

    region = params['region']

    # sequences = ['MT862689']
    region = 'glycoprotein G'
    nucleotide_or_codon = 'nucleotide'
    start_coordinate = ''
    end_coordinate = ''
    reference_sequence = None

    alignment = Alignment(database="RABV_NEW", 
                                sequences=sequences,
                                region=region, 
                                nucleotide_or_codon=nucleotide_or_codon, 
                                start_coordinate=start_coordinate, 
                                end_coordinate=end_coordinate)
    
    data = alignment.get_alignments()
    alignment.build_alignment_fasta_file(data)
    


    with open('my_fasta.fasta', 'r') as file:
        response = HttpResponse(file, content_type='text')
        response['Content-Disposition'] = 'attachment; filename=my_fasta.fasta'
        # os.remove('my_fasta.fasta')

    return response


@api_view(['GET'])
def run_blast_search(request, query):
    tasks_helper = Tasks(database='RABV_NEW')
    queue = django_rq.get_queue('default')  # Get the default queue
    job = queue.enqueue(tasks_helper.runSequenceAlignment)
    
    

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


# OLD ALIGNMENT DOWNLOAD
    # placeholders = ', '.join(['%s'] * len(sequences))  # Creates "%s, %s" dynamically
    # query = f'SELECT s.* FROM sequence s WHERE s.sequence_id IN ({placeholders});'
    


    # with connections[database].cursor() as cursor:
    #     # cursor.execute('SELECT f.ref_seq_name, f.ref_start, f.ref_end, s.* FROM sequence s LEFT JOIN features f on s.alignment_name = f.ref_seq_name WHERE f.product = %s AND f.feature_name="CDS" AND s.sequence_id IN (%s);', [feature_region, sequences])
    #     cursor.execute(query, sequences)
    #     alignments = dictfetchall(cursor)
    #     # print(alignments)
    #     cursor.execute("SELECT ref_start, ref_end FROM features WHERE ref_seq_name='NC_001542' AND feature_name='CDS' AND product=%s", [region])

    #     master_alignment = dictfetchall(cursor)


    # ofile = open("my_fasta.txt", "w")

    
    # ref_start = master_alignment[0]["ref_start"]
    # ref_end = master_alignment[0]["ref_end"]
    # codons = get_kuiken2006_codon_labeling(ref_start, ref_end)
    # codon_start = codons[0]
    # codon_end = codons[1]
    # if (params["coordinates"] == "nucleotide"):

    #     if(params["start"] != ''):
    #         ref_start = int(params["start"])
    #     if(params["end"] != ''):
    #         ref_end = int(params["end"])
    # else:
    #     if(params["start"] != ''):
    #         codon_start = int(params["start"])
    #     if(params["end"] != ''):
    #         codon_end = int(params["end"])
    
    # for i in range(len(alignments)):
    #     if (params["coordinates"] == "nucleotide"):
    #         sub_seq = alignments[i]["alignment"][ref_start:ref_end+1]
    #     else:
    #         sub_seq = alignments[i]["alignment"][ref_start:ref_end+1]
    #         codons = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]
    #         selected_codons = codons[codon_start-1:codon_end]
    #         sub_seq = ''.join(selected_codons)

    #     if set(sub_seq) != {"-"}:
    #         ofile.write(">" + alignments[i]["sequence_id"] + "\n" + sub_seq + "\n")

    # # #do not forget to close it
    # ofile.close()