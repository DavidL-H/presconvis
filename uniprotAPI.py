# Access uniprot queries
import requests, sys, json

WEBSITE_API = "https://rest.uniprot.org/"
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api/"

# Helper function for get response from requested URL
def get_url(url, **kwargs):
    response = requests.get(url, **kwargs);

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit
    
    return response

# Basic search request, searching all UniprotKB
r = get_url(WEBSITE_API+"/uniprotkb/search?query=*")
# The response looks good (code: 200)
print(r)
# Lets have a look at the data:
data = r.json()
n_results = len(data["results"])
print("Number of results: "+ str(n_results))
# There are only 25 results, beacuase search results are paginated
# e.g. this query only return s the first page of results = 25

for (key,value) in r.headers.items():
    print(key+": "+value)

# If you look at the "Link:" key:value pair, you can see the new page in the query, 
# and you can use this to pages iteratively
r.headers["Link"]

