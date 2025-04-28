# Viral Genome Database Resource 
Viral Genome Database Resource (BACKEND)

This codebase is the backend that communitcates with the databases produced by the Viral Genome Toolkit (V-gTK) as well as the Viral Genome Database (V-gDB) web resources. 

It contains the logic for the API that interacts with the database

## Viral Genome Database Structure 
``` line_numbers

V-gDB/
├── api/            # Contains API Endpoints
│   ├── routes/     
│   ├── apps.py
│   └── urls.py     
├── backend/        # Settings for Django project
│   ├── settings.py
│   ├── views.py
│   └── urls.py
├── frontend/       # HTML scripts for V-gDB frontend
│   ├── example_code_scripts/
│   └── static/
├── models/         # Main model logic (including queries)
├── tests/
├── manage.py     
├── .gitignore   
├── README.md     
```

## Local Development

V-gDB is a 

### Requirements:

### Installation:

Install pip (if not installed)

``` bash
# Creating virtual environment

```
``` bash
# Starting virtual environment

```

In your virtual environment:
``` bash
# Install django
python -m pip install Django
```
``` bash
# clone v-gDB app from CVR Github
git clone git@github.com:dana-allen/backend.git
```



### Settings
You need to connect your local database to the vDB-backend app. Inside backend/settings.py, update the following lines with your local database information:

``` python line_numbers
91 DATABASES = {
92    
93    'default': {
94        'ENGINE': 'django.db.backends.mysql', # replace with the type of database
95        'NAME': 'GLUE_TOOLS', # replace with name of database or location of file
96        'USER': 'gluetools', # replace with database user
97        'PASSWORD': 'glue12345', # replace with password
98        'HOST': 'localhost', # replace with host
99        'PORT':'3306', # replace with port
100    }
101 }
```

## Running the project
``` bash
cd vDB-backend # cd into main vDB-backend directory
py manage.py runserver 
```




## Modules Structure 


### API 
API contains all the vDB endpoints available to access the database from a web resource or local script.

1. urls.py

    Points to the routes/urls.py folder. If you want to add more url endpoints that are not located in the routes/urls.py, the file location needs to be added to this file.  

2. apps.py

3. routes

    1. urls.py
        All the url endpoints available. If you are adding more endpoints, each one must be declared in this file. 

        example:
    ``` python
    path('sequences/get_sequences_meta_data/', sequences.get_sequences_meta_data, name='get_sequences_meta_data')

    equivalent to 

    localhost:8000/api/sequences/get_sequences_meta_data/
    ```


### Backend
1. settings.py


``` python line_numbers
61 CORS_ALLOWED_ORIGINS = [
62     "http://localhost:3000",  # Adjust based on your frontend URL
63 ]
```




2. urls.py

### Models

### Tests
