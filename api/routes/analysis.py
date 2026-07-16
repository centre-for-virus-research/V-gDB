from rest_framework.decorators import api_view
from rest_framework.response import Response

from models.helpers import *
from models.helpers import *
from rest_framework import status
from urllib.parse import unquote

from io import StringIO
import uuid
from pathlib import Path
from Bio import SeqIO
import shutil

from models.vgtk_analysis import phylogenetic_clade_assignment_analysis

BASE_DIR = Path(__file__).resolve().parent.parent.parent  # goes up from api/
TASKS_DIR = BASE_DIR / "tasks"

@api_view(['POST'])
def run_phylogenetic_clade_assignment_analysis(request):

    database = request.headers.get('database', 'default')
    file = request.FILES.get("file")
    fasta_text = request.data.get("fasta")

    if file:
        fasta_text = file.read().decode("utf-8")

    if not fasta_text:
        return Response({"error": "No FASTA provided"}, status=400)

    fasta_text = "\n".join(
        line.split()[0] if line.startswith(">") else line
        for line in fasta_text.splitlines()
    )
    # # # Unique job id
    # job_id = "d23a0b10-1335-4668-b351-2c954e7c55df"
    job_id = str(uuid.uuid4())
    job_dir = TASKS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    # # # Folder to hold input FASTA
    inputs_dir = job_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)  # creates tasks/job_id/inputs

    results_dir = job_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)  # creates tasks/job_id/inputs

    # # Wrap fasta_text in StringIO so SeqIO can read it
    fasta_io = StringIO(fasta_text)

    # Parse the input fasta and write each sequence to its own file
    for record in SeqIO.parse(fasta_io, "fasta"):
        seq_id = record.id  # take the FASTA header id
        # Optional: sanitize seq_id to be filename-safe
        safe_id = "".join(c if c.isalnum() or c in "-_." else "_" for c in seq_id)
        seq_file = inputs_dir / f"{safe_id}.fa"

        # Write single sequence FASTA
        SeqIO.write(record, seq_file, "fasta")

    print(f"Split sequences written to {inputs_dir}")
    # try:
    input_path = job_dir / "inputs" /"input.fa"

    with open(input_path, "w") as f:
        f.write(fasta_text)

    try:
        # Run your analysis
        data = phylogenetic_clade_assignment_analysis(database=database, job_id=job_id)

    except ValueError as e:
        print(f"Error: {e}")
        return Response(
            {"error": str(e)},
            status=status.HTTP_400_BAD_REQUEST
        )

    finally:
        # Clean up everything after analysis
        if job_dir.exists() and job_dir.is_dir():
            shutil.rmtree(job_dir)
            print(f"Deleted temporary job folder: {job_dir}")

    # Return the response after successful analysis
    return Response({
        "job_id": job_id,
        "results": data
    })


