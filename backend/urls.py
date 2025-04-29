"""
URL configuration for backend project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from .views import home, installation, schema, check_db_connection, api, mutations, examples, python_tools, about

urlpatterns = [
    path('admin/', admin.site.urls),
    path('api/', include('api.urls')),
    # path('api/', include('api.sequences.urls')),  
    # path('api/', include('api.helpers.urls')),
    # path('api/', include('api.filters.urls')),
    # path('api/', include('api.routes.urls')),

    path('django-rq/', include('django_rq.urls')),
    # path('api/', include('api.urls_hcv'))

    path('', home, name='home'),
    path('installation/', installation, name='installation'),
    path('schema/', schema, name='schema'),
    path('api/', api, name='api'),
    path('mutations/', mutations, name='mutations'),
    path('examples', examples, name='examples'),
    path('python_tools', python_tools, name='python_tools'),
    path('about', about, name='about'),
    path('check-db-connection/', check_db_connection, name='check_db_connection'),
]
