from rest_framework.decorators import api_view
from rest_framework.response import Response


from models.features import Features


@api_view(['GET'])
def get_features(request):

    database = request.headers.get('database', 'default')

    features_helper = Features(database=database)
    features = features_helper.get_features()

    return Response(features)

