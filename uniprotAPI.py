# Accessing the uniprot API with python requests
# Inspired by the following presentation: https://drive.google.com/file/d/1qZXLl19apibjszCXroC1Jg2BGjcychqX/view
# David Lennox-Hvenekilde 28/7/22

import requests, sys, json, time, os
import pandas as pd

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


# Lets do a loop for accessing the entirety of the uniprot query pages ######################################
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

# Turn the response into a dataframe of values of interest #####################################################
def response_to_df(response):
    """
    Turns response into JSON dict, then into a pandas df, using the key values listed in "names".
    """
    all_data = response.json()["results"]
    names = ["scientificName","lineage","entryType","primaryAccession","uniProtkbId","fullName","sequence"]
    newdict = {}
    temp_list1 = []
    temp_list2 = []
    temp_list3 = []
    temp_list4 = []
    temp_list5 = []
    temp_list6 = []
    temp_list7 = []
    for i in range(len(all_data)):
        temp_list1.append(all_data[i]["organism"][names[0]])
        temp_list2.append(", ".join(all_data[1]["organism"][names[1]]))
        temp_list3.append(all_data[i][names[2]])
        temp_list4.append(all_data[i][names[3]])
        temp_list5.append(all_data[i][names[4]])
        temp_list6.append(all_data[i]["proteinDescription"]["recommendedName"][names[5]]["value"])
        temp_list7.append(all_data[i][names[6]]["value"])
    list_of_lists = [temp_list1,temp_list2,temp_list3,temp_list4,temp_list5,temp_list6,temp_list7]
    #length = [len(temp_list1), len(temp_list2),len(temp_list3), len(temp_list4), len(temp_list5), len(temp_list6), len(temp_list7)]

    i=0
    for n in names:
        newdict[n] = list_of_lists[i]
        i += 1

    query_df_temp = pd.DataFrame(newdict)
    return query_df_temp

# Grab all query pages as dataframe ###################################################################
def all_df_query(search_string, save_csv = False):
    """
    Takes a Uniprot query string as input:
    e.g. "uniref_cluster_90:UniRef90_Q8X825"
    And returns a pandas df with key values (defined by response_to_df())
    of all the entries in the query
    """
    # Start the loop to extract ID and sequence from the JSON query response
    start = time.time()
    r = uniprot_query(search_string)

    # While there is still a next page:
    inf_loop = True
    n_page = 1
    while inf_loop:
        # Extract key values from json to dataframe
        try:
            df_temp = response_to_df(r)
        except:
            break

        # Try to append temporary dataframe to master dataframe
        try: 
            frames = [query_df,df_temp]
            query_df = pd.concat(frames)
        except:
            query_df = df_temp

        # Try to grab the next page of the query results, if no next page there, break.
        try:
            next_page = r.headers["Link"]
        except:
            break

        # Extract next page url, and use for querying next page of search
        next_page = next_page[next_page.index("<")+1:next_page.index(">")]
        r = get_url(next_page)
        n_page += 1

    end = time.time()
    print("dataframe fetched in  "+ str(round(end - start, 3)) + " seconds")

    if save_csv:
        query_df.to_csv(search_string[search_string.index(":")+1::]+".csv")
    else:
        return query_df