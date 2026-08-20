from rest_framework.decorators import api_view
from rest_framework.response import Response


@api_view(['GET'])
def get_version(request):

    return Response("V-gDB has successfully been installed!")
