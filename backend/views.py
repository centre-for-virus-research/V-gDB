from django.db import connections
from django.db.utils import OperationalError
from django.shortcuts import render
from django.urls import get_resolver, get_urlconf
import importlib
from django.conf import settings
from django.http import JsonResponse
from collections import defaultdict
import json
def home(request):
    """Renders a styled home page with available API endpoints from api/routes/urls.py."""
    
    # Import your specific URL configuration
    api_urls_module = importlib.import_module("api.urls")
    url_patterns = getattr(api_urls_module, "urlpatterns", [])

    # Extract endpoint paths
    endpoints = []
    for pattern in url_patterns:
        if hasattr(pattern, 'pattern') and hasattr(pattern.pattern, '_route'):
            endpoints.append("/api/"+pattern.pattern._route)

    # try:
    #     # Try to connect to the database and check if the connection is available
    #     connections['default'].ensure_connection()
    #     db_connected = True
    #     # If the connection is successful, you can return a success message or render a template
    # except OperationalError:
    #     # If there's an error (no database connection), handle it
    #     db_connected = False

    db_connected = True
    db_name = ''
    try:
        # Attempt to fetch a connection to the database
        connection = connections['default']
        db_name = connection.settings_dict['NAME']
        # This will trigger the connection to the database and check if it's live
        with connection.cursor():
            pass
    except OperationalError:
        db_connected = False
        db_name = "Unknown (Connection Failed)"

    db_name = settings.DATABASES['default']['NAME']
    return render(request, 'index.html', {'endpoints': endpoints, 
                                          'debug': settings.DEBUG, 
                                          'db_connected':db_connected, 
                                          "db_name":db_name})


def installation(request):
    """Renders installation guide."""


    return render(request, 'installation.html')

def schema(request):
    """Renders database schema."""


    return render(request, 'schema.html')

def python_tools(request):
    """Renders database schema."""

    return render(request, 'python_tools.html')

def about(request):
    """Renders database schema."""

    return render(request, 'about.html')

def mutations(request):
    """Renders mutations schema."""

    return render(request, 'mutations.html')

def examples(request):
    """Renders examples schema."""

    return render(request, 'examples.html')

def gui(request):
    """Renders gui schema."""

    return render(request, 'gui.html')

def api(request):
    """Renders API schema grouped by top-level type."""

    api_urls_module = importlib.import_module("api.urls")
    url_patterns = getattr(api_urls_module, "urlpatterns", [])

    # Group endpoints by category
    endpoints = defaultdict(list)

    for pattern in url_patterns:
        if hasattr(pattern, 'pattern') and hasattr(pattern.pattern, '_route'):
            route = pattern.pattern._route
            full_path = "/api/" + route
            category = route.split('/')[0] if '/' in route else 'uncategorized'
            endpoints[category].append(full_path)

    with open("/Users/danaallen/CVR/gdb/web-resources/V-gDB/frontend/static/javascript/openapi_new.json", 'r') as f:
        data = json.load(f)

    return render(request, 'api.html', {'endpoints': data})

def check_db_connection(request):
    try:
        # Try to connect to the database
        connections['default'].cursor()
        return JsonResponse({'status': 'success'})
    except OperationalError:
        # If there's an error, the database is not connected
        return JsonResponse({'status': 'failure'})