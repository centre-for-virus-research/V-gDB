from rest_framework.decorators import api_view
from rest_framework.response import Response
from models.mutations import Mutations


@api_view(['GET'])
def get_host_adaptations(request):

    database = request.headers.get('database', 'default')
    params = dict(request.GET.items())

    final_params = {}
    master_accession = ""
    if params:
        for key, value in params.items():
            params[key] = value.split(',') if ',' in value else value
        if "phylum" in params:
            final_params = {"phylum":params["phylum"]}
        if "class" in params:
            final_params = {"class":params["class"]}
        if "family" in params:
            final_params = {"family":params["family"]}
        if "genus" in params:
            final_params = {"genus":params["genus"]}
        if "species" in params:
            final_params = {"species":params["species"]}
        if "master_accession" in params:
            master_accession = params["master_accession"]
            del params["master_accession"]

    mutations = Mutations(database=database)

    data = None
    if database == "RABV" and params:
        # data = mutations.get_adaptive_mutations_chart_RABV(final_params)
        data = mutations.get_host_adaptations(final_params, master_accession)
    # else :
    #     data = mutations.get_adaptive_mutations_chart(segment)

    return Response(data)
