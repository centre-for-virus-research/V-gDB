from django.shortcuts import render
from rest_framework.decorators import api_view
from rest_framework.response import Response

from models.helper_functions import *

import django_rq

from models.tasks import Tasks
from models.alignment import Alignment
from models.features import Features
from django.http import HttpResponse
from Bio import SeqIO 

@api_view(['GET'])
def get_blast_results(request, file_path):

    with open("/Users/dana/CVR/backend/jobs/"+file_path+"/blast_results.tsv", 'r') as file:
        response = HttpResponse(file, content_type='tsv')
        response['Content-Disposition'] = 'attachment; filename=blast_results.tsv'
    return response

@api_view(['GET'])
def get_nextAlign_results(request, file_path):
    print("WE ARE STARTING THIS")
    sequences = []
    file_path = "/Users/dana/CVR/backend/jobs/"+file_path+"/Nextalign/MT862689/MT862689.aligned.fasta"
    # Read FASTA file and store sequences
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append({"query": record.id, "seq": str(record.seq)})

    if not sequences:
        return None  # Return None if the file is empty

    # Extract primary accession (first sequence)
    primary_accession = sequences[0]["query"]
    primary_seq = sequences[0]["seq"]

    # Extract aligned sequences (everything except the first one)
    aligned_sequences = sequences[1:]

    features_helper = Features(database='RABV_NEW')

    features = features_helper.get_feature(primary_accession) 


    # Construct final JSON output
    result = [{
        "primary_accession": primary_accession,
        "seq": primary_seq,
        "alignedSeq": aligned_sequences,
        "features": features
    }]
    print("HIHIHI")
    print(result)
    return Response(result)


@api_view(['GET'])
def run_sequence_alignment(request):
    params = dict(request.GET.items())
    query = params["query"]
    alignment_helper = Alignment(database='RABV_NEW', sequences=query)

    try:
        queue = django_rq.get_queue('default')  # Get the default queue
        job = queue.enqueue(alignment_helper.runSequenceAlignment)
        return Response(job.id)
    except ConnectionError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)



@api_view(['GET'])
def get_job_logs(request, job_id):
    queue = django_rq.get_queue('default')
    job = queue.fetch_job(job_id)
    if job and job.meta.get('status'):
        return Response(job.meta)
    return Response({"status": []})

@api_view(['GET'])
def get_job_result(request, job_id):
    queue = django_rq.get_queue('default')  # Replace 'default' with your queue name if needed
    job = queue.fetch_job(job_id)

    if job is None:
        return Response({"error": "Job not found", "status":404})

    if job.is_failed:
        return Response({"status": "failed", "error": str(job.exc_info)})

    return Response({"status": str(job.get_status()), "result": str(job.result)})

