from django.db import connections
import csv
from models.helpers import *
from collections import Counter
from collections import defaultdict

class Strains:
    """
    A class to handle operations related to sequence metadata, alignments, 
    and filtered retrieval from a database in a Django project.
    """

    def __init__(self, database, filters=None):
        """
        Initialize the Sequences instance.

        Args:
            database (str): The name of the database alias to connect to.
            primary_accession (str, optional): Identifier for a specific sequence. Defaults to None.
            filters (dict, optional): Dictionary of filters for querying metadata. Defaults to None.
        """
        self.database = database  
        self.filters = filters

    def get_strains(self):
        # if not self.filters:
        #     with connections[self.database].cursor() as cursor:
        #         cursor.execute(' SELECT isolate, country, host, segment, collection_date, primary_accession FROM meta_data ORDER BY isolate, segment;')
        #         rows = cursor.fetchall()
        # # Structure data by isolate
        # data_by_isolate = defaultdict(lambda: {"host": None, "segments": {}})

        # for isolate, country, host, segment, collection_date, accession in rows:
        #     if data_by_isolate[isolate]["host"] is None:
        #         data_by_isolate[isolate]["host"] = host
        #         data_by_isolate[isolate]["country"] = country
        #         data_by_isolate[isolate]["collection_date"] = collection_date
        #     data_by_isolate[isolate]["segments"][segment] = accession

        # # Convert to a nice list/dict if needed
        # result = []
        # for isolate, data in data_by_isolate.items():
        #     entry = {
        #         "isolate": isolate,
        #         "country": data["country"],
        #         "collection_date": data["collection_date"],
        #         "host": data["host"],
        #         "segments": data["segments"]  # dict like {1: 'ACC001', 2: 'ACC002', ...}
        #     }
        #     result.append(entry)

        query = f"""
                    SELECT i.*, md.host, md.collection_year, md.country 
                    FROM isolates i
                    JOIN meta_data md ON md.strain = i.strain
                """
        with connections[self.database].cursor() as cursor:
                cursor.execute(query)
                result = dictfetchall(cursor)

        return result

    def get_strain(self, isolate):
        if isolate:
            with connections[self.database].cursor() as cursor:
                cursor.execute('SELECT * FROM meta_data where isolate=%s ORDER BY segment', [isolate])
                result = dictfetchall(cursor)
                # print(result)
                for l in result:
                    print(l["primary_accession"])
                    print("here0")
                    cursor.execute("SELECT * FROM sequence_alignment WHERE sequence_id=%s", [l["primary_accession"]])
                    alignment = dictfetchall(cursor)
                    print("here4")
                    if alignment:
                        l["alignment"] = alignment[0]
                        print("here")
                        # Add insertions
                        cursor.execute("SELECT * FROM insertions WHERE accession = %s", [l["primary_accession"]])
                        l["alignment"]["insertions"] = dictfetchall(cursor)
                        print("here2")
                        # Add features
                        print(l["alignment"]["alignment_name"])
                        cursor.execute("SELECT * FROM features WHERE accession=%s and reference_accession = %s ORDER BY cds_start", [l["primary_accession"], l["alignment"]["alignment_name"]])
                        l["alignment"]["features"] = dictfetchall(cursor)
                        print("here3")

                        # Add reference sequence
                        cursor.execute("SELECT alignment FROM sequence_alignment WHERE sequence_id = %s", [l["alignment"]["alignment_name"]])
                        l["alignment"]["ref_seq"] = dictfetchall(cursor)[0]["alignment"]

                # print(result)
                
        
        
        return result
