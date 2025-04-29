from rest_framework.decorators import api_view
from rest_framework.response import Response


from models.features import Features


@api_view(['GET'])
def get_features(request):

    database = request.headers.get('database', 'default')

    features = Features(database=database)
    features = features.get_features()

    return Response(features)

