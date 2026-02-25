import requests

# url = "http://gdb-dev.cvr.gla.ac.uk/api/sequence/MT862689"
url = "http://localhost:8000/api/sequence/MT862689"

headers = {
    "Content-Type": "application/json",
    'database': "FLUV"
}

response = requests.get(url, headers=headers)

print("Status code:", response.status_code)
print("Response:", response.json())