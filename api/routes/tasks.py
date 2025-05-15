from rest_framework.decorators import api_view
from rest_framework.response import Response
import django_rq
from django.http import HttpResponse, Http404
from Bio import SeqIO
import os

from models.helpers import *
from models.tasks import Tasks
from models.features import Features

# ----------------------------------------------------------------------
# View to enqueue a sequence alignment job
#     # redis-server IN THE COMMAND LINE
#     # python manage.py rqworker default MUST RUN THIS TO GET IT TO WORK
# ----------------------------------------------------------------------
@api_view(['GET'])
def run_sequence_alignment(request):
    """
    Enqueues a sequence alignment task using RQ (Redis Queue).
    Requires a 'query' parameter and optionally a 'database' header.
    Returns:
        Job ID (str): The ID of the queued job.
    """
    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())
    query = params.get("query")

    tasks = Tasks(database=database, query=query)

    try:
        queue = django_rq.get_queue('default')
        job = queue.enqueue(tasks.run_sequence_alignment)  # Submit task to queue
        return Response(job.id)
    except ConnectionError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
    

# ----------------------------------------------------------------------
# View to get sequence alignment results in FASTA format
# ----------------------------------------------------------------------
@api_view(['GET'])
def get_alignment_results(request, job_id):
    """
    Parses aligned FASTA sequences and metadata for a given job ID.
    Args:
        job_id (str): ID of the job to retrieve alignment results for.
    Returns:
        Response: JSON containing primary sequence, aligned sequences, and features.
    """

    sequences = []
    file_path = os.path.join('jobs', job_id, 'analysis','query_ref_alignment')

    if not os.path.exists(file_path):
        raise Http404("Alignment file not found.")
    for f in os.listdir(file_path):

        with open(os.path.join(file_path, f), "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):

                sequences.append({"query": record.id, "seq": str(record.seq)})

    if not sequences:
        return Response([])  # Empty response if no sequences found

    primary_accession = sequences[0]["query"]
    primary_seq = sequences[0]["seq"]
    aligned_sequences = sequences[1:]

    features_helper = Features(database='RABV')
    features = features_helper.get_feature(primary_accession)

    result = [{
        "primary_accession": primary_accession,
        "seq": primary_seq,
        "alignedSeq": aligned_sequences,
        "features": features
    }]

    return Response(result)

# ----------------------------------------------------------------------
# View to download BLAST results as a TSV file
# ----------------------------------------------------------------------
@api_view(['GET'])
def get_blast_results(request, job_id):
    """
    Downloads the BLAST results file for a given job ID.
    Args:
        job_id (str): ID of the job to retrieve results for.
    Returns:
        HttpResponse: File response with TSV data.
    """
    file_path = os.path.join('jobs', job_id, 'analysis', 'query_tophits_uniq.tsv')

    if not os.path.exists(file_path):
        raise Http404("BLAST results file not found.")

    with open(file_path, 'r') as file:
        response = HttpResponse(file, content_type='text/tab-separated-values')
        response['Content-Disposition'] = 'attachment; filename=blast_results.tsv'
        return response



# ----------------------------------------------------------------------
# View to fetch job logs or metadata for a queued task
# ----------------------------------------------------------------------
@api_view(['GET'])
def get_job_logs(request, job_id):
    """
    Retrieves metadata (including status) of a job from the RQ queue.
    Args:
        job_id (str): Job ID to fetch logs for.
    Returns:
        Response: Job metadata dictionary or empty status if not found.
    """
    queue = django_rq.get_queue('default')
    job = queue.fetch_job(job_id)

    if job and job.meta.get('status'):
        return Response(job.meta)

    return Response({"status": []})
