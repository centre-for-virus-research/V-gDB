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

If you wish to exapnd upon V-gDB functionality or API, you can download V-gDB and install locally.

### 1. Installation:

``` bash
# clone v-gDB app from CVR Github
git clone git@github.com:centre-for-virus-research/V-gDB.git 
```
``` bash
cd V-gDB 
```
``` bash
# Creating virtual environment
python3 -m venv V_gDB_env 
```
``` bash
# Starting virtual environment
source V_gDB_env/bin/activate 
```
``` bash
# Installing V-gDB requirements
pip3 install -r requirements.txt 
```
``` bash
python3 manage.py migrate 
```

### 2. Connecting Database:
You need to connect your local database to the vDB-backend app. You can download one of V-gDB's pre-compiled databases, or generate your own using the Viral Genome Toolkit.

Once you have the .db file, inside backend/settings.py, update the following lines with your local database information:

``` python line_numbers
DATABASES = {
    'default': {
    'ENGINE': 'django.db.backends.sqlite3',
    'NAME': 'db/rabv.db', #Update with local path to local database .db file
    'HOST': 'localhost', 
    'PORT':'3306', 
    }
}
```

### 3. Starting V-gDB:

If the V-gDB app is installed and the database information has been updated, run this in your V-gDB directory:
``` bash
# Starting V-gDB server
python3 manage.py runserver
```
Navigate to https://localhost:8000. If you see the following screen, you have successfully installed V-gDB!

![success](frontend/static/images/success_install.png?raw=true)

## Modules Structure 


### API 
API contains all the vDB endpoints available to access the database from a web resource or local script.

1. urls.py

    Points to the routes/urls.py folder. If you want to add more url endpoints that are not located in the routes/urls.py, the file location needs to be added to this file.  

2. apps.py


### Models

### Tests
