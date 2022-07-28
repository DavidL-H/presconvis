# Access uniprot queries
import requests, sys, json, time

# rest end-points
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

##################################################################
# Lets build a function:
# the query is what comes after https://www.uniprot.org/uniprotkb?query=
# When you are searching manually on the website

search_string = "uniref_cluster_90:UniRef90_Q8X825"
def uniprot_query(search_string):
    """
    This function returns the json response of a advanced search query.
    This search query is what comes after 'https://www.uniprot.org/uniprotkb?query='
    when doing queries on the uniprot website.
    e.g.
    uniprot_query('uniref_cluster_90:UniRef90_Q8X825')
    """
    r = get_url(WEBSITE_API+"/uniprotkb/search?query="+search_string)
    return r

def uniprot_json_to_fasta_print(json_response):
    """
    This takes a json_response from get_url() on the uniprotkb API, and returns all the 
    protein IDs and proteins sequences printed in the console
    """
    data = json_response.json()
    for i in range(len(data["results"])):
        print(">"+data["results"][i]["primaryAccession"])
        print(data["results"][i]["sequence"]["value"])


def write_fasta_from_json(json_response, file_name):
    """
    This takes a json_response from get_url() on the uniprotkb API, and a search string (file_name) returns all the 
    protein IDs and proteins sequences and saves them as a .fasta file
    """
    data = json_response.json()
    for i in range(len(data["results"])):
        with open(file_name[file_name.index(":")+1::]+".fasta", "a") as fasta_file:
            fasta_file.writelines(">"+data["results"][i]["primaryAccession"]+"\n")
            fasta_file.writelines(data["results"][i]["sequence"]["value"]+"\n")


# Lets do a loop for accessing the entirety of the uniprot query pages
def all_fasta_query(search_string):
    """
    This function loops through all the pages of a uniprot query, saving all entries to a .fasta file.
    e.g.
    search_string = "uniref_cluster_90:UniRef90_Q8X825"
    all_fasta_query(search_string)
    """
    start = time.time()
    r = uniprot_query(search_string)

    # While there is still a next page:
    inf_loop = True
    n_page = 1
    while inf_loop:
        # Write the response to fasta
        write_fasta_from_json(r, search_string)
        try:
            next_page = r.headers["Link"]
        except:
            #print("last page")
            break
        next_page = next_page[next_page.index("<")+1:next_page.index(">")]
        #print(next_page)
        r = get_url(next_page)
        n_page += 1

    end = time.time()
    print("fasta file fetched in "+ str(round(end - start, 3)) + " seconds")

# Example for downloading all sequences of uniref_90 cluster Q8X825 (BioB)
search_string = "uniref_cluster_90:UniRef90_Q8X825"
all_fasta_query(search_string)

