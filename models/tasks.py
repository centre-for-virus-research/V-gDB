import logging
import rq
import os
from os.path import join
from pathlib import Path


logger = logging.getLogger(__name__)

# ALL OF THESE USE redis server to work
class Tasks:
    def __init__(self, database, query):
        self.database = database  
        self.query = query

    @staticmethod
    def path_to_basename(file_path):
        path = os.path.basename(file_path)
        return path.split('.')[0]


    def process_genbank_submission(self):
        db_path = Path(join("analysis","db.fa"))
        job = "test"
        os.makedirs(join("jobs",job), exist_ok=True)

    # OFFICEAL SEQUENCE ALIGNMNET TASK CODE
    def run_sequence_alignment(self):
        # Add in sequence alignment logic here
        logger.info(f"Starting sequence alignment")
        job = rq.get_current_job()
        job.meta['status'] = "done"
        job.save_meta()

        return job.meta


