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


from models.sequences import Sequences
from models.mutations import Mutations

@api_view(['GET'])
def get_sequences_meta_data(request):

    database = request.headers.get('database', 'default')

    sequences = Sequences(database=database)

    try:
        data = sequences.get_sequences_meta_data()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
        
    return Response(data)

@api_view(['GET'])
def get_sequence_meta_data(request, primary_accession):

    database = request.headers.get('database', 'default')

    sequences_helper = Sequences(database=database, primary_accession=primary_accession)

    try:
        data = sequences_helper.get_sequence_meta_data()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)
    
    # primary_accession = 'MT862689'
    # mutations_helper = Mutations(database="RABV_NEW",
                                # reference_sequence='EF437215'
                                # )
    # f = open("/Users/dana/CVR/backend/models/ef437215.txt", "r")
    # test_sequence = f.read()
    # data["mutations"] = mutations_helper.get_sequence_mutations(primary_accession, test_sequence)

    return Response(data)

@api_view(['GET'])
def get_sequences_meta_data_by_filters(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())

    for key, value in params.items():
        params[key] = value.split(',') if ',' in value else value

    sequences_helper = Sequences(database=database, filters=params)

    try:
        data = sequences_helper.get_sequences_meta_data_by_filters()
    except ValueError as e:
        print(f"Error: {e}")
        return HttpResponse(e, status=404)

    return Response(data)

@api_view(['GET'])
def download_sequences_meta_data(request):

    database = request.headers.get('database', 'RABV_NEW')
    params = dict(request.GET.items())

    sequences = Sequences(database=database, filters=params)
    if params:
        data = sequences.get_sequences_meta_data_by_filters()
    else:
        data = sequences.download_sequences_meta_data()
    sequences.build_meta_data_csv_file(data)

    with open('tmp_file.csv', 'r') as file:
        response = HttpResponse(file, content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename=tmp_file.csv'
        os.remove('tmp_file.csv')

        return response
    
# THIS IS BEING USED IN THE MUTATION GUI
# TODO: REMOVE this and use callback version to get host
@api_view(['GET'])
def get_host_species(request):
    database = request.headers.get('database', 'default')
    
    with connections[database].cursor() as cursor:
        cursor.execute('SELECT DISTINCT(host) FROM meta_data WHERE host IS NOT NULL;')
        result = dictfetchall(cursor)
    return Response(result)



@api_view(['GET'])
def get_filtered_sequences(request):
    database = request.headers.get('database', 'default')
    where_clauses = []
    params = []
    country_clauses = []
    country_params = []
    
    params = dict(request.GET.items())
    for key, value in params.items():
        params[key] = value.split(',') if ',' in value else value

    for key, value in params.items():
        new_value = value
        # new_value = value.split(',') if ',' in value else value
        # new_value = value

        if key == 'gb_length_lower':
            where_clauses.append(f"length >= %s ")
            params.append(new_value)
        elif key == 'gb_length_upper':
            where_clauses.append(f"length <= %s ")
            params.append(new_value)

        elif key == 'gb_update_date_lower':
            where_clauses.append(f"YEAR(update_date) >= %s ")
            params.append(new_value)
        elif key == 'gb_update_date_upper':
            where_clauses.append(f"YEAR(update_date) <= %s ")
            params.append(new_value)

        elif key == 'creation_year_lower':
            where_clauses.append(f"YEAR(create_date) >= %s ")
            params.append(new_value)
        elif key == 'creation_year_upper':
            where_clauses.append(f"YEAR(create_date) <= %s ")
            params.append(new_value)

        elif key == 'collection_year_lower':
            where_clauses.append(f'STRFTIME("%Y", collection_date) >= %s ')
            params.append(new_value)
        elif key == 'collection_year_upper':
            where_clauses.append(f"YEAR(collection_date) <= %s ")
            params.append(new_value)

        elif key in ['primary_accession', 'host', 'isolate', 'pubmed_id']:
            new_value = [new_value] if isinstance(new_value, str) else new_value
            where_clauses.append(f"{key} IN ({','.join(["%s"] * len(new_value))})")
            params.extend(new_value)

        elif key in ['major_clade', 'minor_clade']:
            where_clauses.append(f"{key} IN %s")
            


        elif key in ['m49_region_id', 'm49_sub_region_id', 'm49_code', 'development_status']:
            country_clauses.append(f"{key} IN %s")
            country_params.append(new_value)

        elif key in ['is_ldc', 'is_lldc', 'is_sids']:
            country_clauses.append(f"{key} = %s")
            country_params.append(True if value == 'true' else False)


    # Country-specific query
    if country_clauses:
        
        country_query = f"SELECT DISTINCT(id) FROM m49_country WHERE {' AND '.join(country_clauses)}"
        with connections[database].cursor() as cursor:
            cursor.execute(country_query, country_params)
            country_ids = [row[0] for row in cursor.fetchall()]
        if country_ids:
            where_clauses.append("m49_country_id IN %s")
            params.append(country_ids)

    # Main query
    where_str = ' AND '.join(where_clauses)
    query = f"SELECT * FROM meta_data \
                WHERE {where_str} ORDER BY create_date DESC;" if where_clauses else "SELECT * FROM meta_data"
    with connections[database].cursor() as cursor:
        cursor.execute(query, params)
        result = dictfetchall(cursor)

    print(result)

        
    return Response(result)

    


@api_view(['GET'])
def advanced_filter(request, query):
    database = request.headers.get('database', 'RABV_NEW')
    print(query)
    query="SELECT * FROM meta_data WHERE"+query

    with connections[database].cursor() as cursor:
        cursor.execute(query)
        result = dictfetchall(cursor)
    return Response(result)



def compute_match_percentage(query, text):
    max_len = max(len(query), len(text))
    if max_len == 0:
        return 100  # Both are empty
    # return round((1 - levenshtein_distance(query, text) / max_len) * 100, 2)


def get_kmers(seq, k=5):
    """Convert sequence into a set of k-mers."""
    return set(seq[i:i+k] for i in range(len(seq) - k + 1))

def jaccard_similarity(seq1, seq2, k=5):
    """Compute Jaccard similarity using k-mers."""
    kmers1 = get_kmers(seq1, k)
    kmers2 = get_kmers(seq2, k)

    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)

    return round((intersection / union) * 100, 2) if union != 0 else 0


@api_view(['GET'])
def search_sequence(request, sequence):

    database = request.headers.get('database', 'RABV_NEW')
    database = 'RABV_NEW'
    with connections[database].cursor() as cursor:
        cursor.execute('SELECT * from sequence')
        result = dictfetchall(cursor)
    search_query = sequence.lower().replace(' ', '')

    matches = []

    for item in result:
        alignment = item['alignment'].replace('-','').replace(' ', '').lower()
        jaccard = jaccard_similarity(alignment, search_query)

        if(jaccard > 95):

            with connections[database].cursor() as cursor:
                cursor.execute('SELECT * from meta_data where primary_accession=%s', [item["sequence_id"]])
                sequence = dictfetchall(cursor)
            item["jaccard_score"] = jaccard
            item["sequence"]=sequence[0]
            matches.append(item)
    
    sorted_matches = sorted(matches, key=lambda x: x["jaccard_score"], reverse=True)


    return Response(sorted_matches)






#alignment:
# - remove gaps in triplets 

# - select a region
# - else select nucleotide range\
# - else select codon range- based off of NC_

# - grab protein region 
#     - if only gaps, ignore it



# @api_view(['GET'])
# def get_mutation_regions_and_codons(request):

#     database = request.headers.get('database', 'RABV_NEW')
#     database = 'RABV_NEW'
#     query = 'SELECT product, ref_start, ref_end FROM features where feature_name="CDS" and ref_seq_name="NC_001542"'
#     with connections[database].cursor() as cursor:
        
#         cursor.execute(query)
#         result = dictfetchall(cursor)
#     for i in result:
#         print("HERERE")
#         print(i)
#         codons = get_kuiken2006_codon_labeling(i["ref_start"], i["ref_end"])
#         i["codon_start"] = codons[0]
#         i["codon_end"] = codons[1]


#     return Response(result)


# @api_view(['GET'])
# def get_mutations(request):

#     database = request.headers.get('database', 'RABV_NEW')
#     alignments = []
#     database = 'RABV_NEW'

#     params = dict(request.GET.items())

    
#     print(params)
#     hosts = params["host"].split(',')

#     placeholders = ', '.join(['%s'] * len(hosts))  
#     query = f'SELECT primary_accession FROM meta_data WHERE host IN ({placeholders})'
#     query = f'SELECT s.*, m.* FROM sequence s LEFT JOIN meta_data m on s.sequence_id = m.primary_accession WHERE m.host IN ({placeholders});'


#     with connections[database].cursor() as cursor:

#         cursor.execute(query, hosts)
#         alignments = dictfetchall(cursor)

#         cursor.execute("SELECT ref_start, ref_end FROM features WHERE ref_seq_name='NC_001542' AND feature_name='CDS' AND product=%s", [params["region"]])

#         master_alignment = dictfetchall(cursor)

#     ref_start = master_alignment[0]["ref_start"]
#     ref_end = master_alignment[0]["ref_end"]
    
#     codon_chosen = params["codon"].split(',')
#     results = []
#     for i in range(len(alignments)):
#         alignments[i]["mutations"] = {}
#         sub_seq = alignments[i]["alignment"][ref_start:ref_end+1]
#         codons = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]

#         for codon in codon_chosen:
#             selected_codons = codons[int(codon)]
#             sub_seq = ''.join(selected_codons)
#             # alignments[i]["mutations"][str(codon)] = {}
#             if "-" not in sub_seq:
#                 if len(sub_seq)%3 == 0:
#                     protein = translate(sub_seq)
#                     alignments[i]["mutations"][str(codon)] = protein
#             else:
#                 if sub_seq != "---":
#                     alignments[i]["mutations"][str(codon)] = "X"


#     return Response(alignments)

# def translate(codon): 
       
#     table = { 
#         'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
#         'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
#         'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
#         'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
#         'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
#         'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
#         'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
#         'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
#         'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
#         'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
#         'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
#         'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
#         'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
#         'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
#         'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
#         'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
#     } 
#     protein ="" 
#     protein = table.get(codon, "-")
#     return protein 