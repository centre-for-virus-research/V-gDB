from django.db import connections
import time
import logging
import django_rq
import django
from itertools import groupby
import rq

import os
import sys

from os.path import join

from pathlib import Path


logger = logging.getLogger(__name__)

# ALL OF THESE USE redis server to work
class Tasks:
    def __init__(self, database):
        self.database = database  

    @staticmethod
    def path_to_basename(file_path):
        path = os.path.basename(file_path)
        return path.split('.')[0]


    def process_genbank_submission(self):
        db_path = Path(join("analysis","db.fa"))
        job = "test"
        os.makedirs(join("jobs",job), exist_ok=True)

    # OFFICEAL SEQUENCE ALIGNMNET TASK CODE
    def run_sequence_alignment(self, query):
        # Add in sequence alignment logic here
        logger.info(f"Starting sequence alignment")
        print(query)
        job = rq.get_current_job()
        job.meta['status'] = "running"

        return job.meta


