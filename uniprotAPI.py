# Accessing the uniprot API with python requests
# Inspired by the following presentation: https://drive.google.com/file/d/1qZXLl19apibjszCXroC1Jg2BGjcychqX/view
# David Lennox-Hvenekilde 28/7/22

import requests, sys, json, time, os

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
    #print(r.headers)
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
def all_fasta_query(search_string, force = True):
    """
    This function loops through all the pages of a uniprot query, saving all entries to a .fasta file.
    Input:
    - search_string: a uniprot advanced query suffix (e.g. "uniref_cluster_90:UniRef90_Q8X825")
    - force: Overwrite the existing fasta file of the same name, if present?
    Output:
    - fasta file saved to current wd, 
    - retuns fasta file path to be saved to variable (e.g. "./UniRef90_Q8X825.fasta")
    """
    # Check is the fasta file already exists in the directory
    file_path = search_string[search_string.index(":")+1::]+".fasta"
    if os.path.isfile(file_path):
        if force:
            os.remove(file_path)
            print("Existing fasta file deleted, and remade")
        elif not force:
            print("Fasta already exists")
            return file_path

    # Start the loop to extract ID and sequence from the JSON query response
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
    return file_path


