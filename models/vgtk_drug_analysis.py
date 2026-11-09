from django.db import connections
from pathlib import Path
from .helpers import *
from Bio.Seq import Seq


BASE_DIR = Path(__file__).resolve().parent.parent # goes up from api/
TASKS_DIR = BASE_DIR / "tasks"
RESOURCES_DIR = BASE_DIR / "resources"

FEATURES_MAP = {'core protein':'Core',
                    'envelope protein E1': 'E1',
                    'envelope protein E2':'E2',
                    'RNA-dependent RNA polymerase NS5B':'NS5B',
                    'nonstructural protein NS5A':'NS5A',
                    'nonstructural protein NS4B':'NS4B',
                    'nonstructural protein NS4A':'NS4A',
                    'protease/helicase protein NS3':'NS3',
                    'nonstructural protein NS2':'NS2',
                    'protein p7':'p7'
                  }

def translate_sequence(sub_seq):
    
    clean_seq = sub_seq.replace("-", "N")
    protein = str(Seq(clean_seq).translate(to_stop=False))

    return protein


def drug_analysis(alignment, features):


    known_mutations = []
    present_mutations = []
    final_mutations = []
    drug_regimen_results = []

    for feature in features:
        clean_feature = FEATURES_MAP[feature['product']]
        if clean_feature not in ["NS3", "NS5A", "NS5B"]:
            continue
        cds_start = int(feature["cds_start"])
        cds_end = int(feature["cds_end"])
        
        print(clean_feature)
        sub_seq = alignment[cds_start-1: cds_end-1]
        protein_seq = translate_sequence(sub_seq)
        # print("PROTEIN ", protein)
        with connections["HCV"].cursor() as cursor:
            cursor.execute("""SELECT distinct (signature_id) 
                                FROM mutation_catalog
                                WHERE protein_name = %s
                            """, [clean_feature])

            mutations = dictfetchall(cursor)
            # print(mutations)

        for mutation in mutations:
            signature_id = mutation["signature_id"]

            protein, signature = signature_id.split(':')

            single_signature = signature.split('+')
            combination = True

            for ss in single_signature:
                if(ss in known_mutations):
                    continue

                aa = ss[-1]
                if "de" in ss[:-1]:
                    combination = False
                    continue
                    
                aa_idx = int(ss[:-1])

                if aa_idx < 1 or aa_idx > len(protein_seq):
                    continue

                query_aa = protein_seq[aa_idx-1]
                if (query_aa == aa):
                    known_mutations.append(ss)
                else:
                    combination = False

            if combination:
                present_mutations.append(signature_id)
                
        
        if len(present_mutations) > 0:
            
            placeholders = ', '.join(['%s'] * len(present_mutations))
            query = f""" SELECT * FROM mutation_catalog WHERE signature_id IN ({placeholders})
                    """
            with connections["HCV"].cursor() as cursor:
                final_mutations.extend(fetch_all(cursor, query, present_mutations))
                
    return {"mutations":final_mutations}
