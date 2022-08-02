from uniprotAPI import all_fasta_query, uniprot_query, all_df_query
from clustalomegaAPI import post_alignment_request, clustalo_alignment, root
#import clustalomegaAPI as clust
import pandas as pd
import time

# Testing out UNIPROT functions
##################################################################

# Test out fetching of uniprot protein sequences from custom function:
#search_string = "uniref_cluster_90:UniRef90_Q8X825"
search_string = "uniref_cluster_50:UniRef50_Q8A7T2"

# Print query metadata
response = uniprot_query(search_string)
print(response.headers)

# Write query to fasta file
fasta_file = all_fasta_query(search_string, force=False)
print(fasta_file)

# Extract key data to dataframe:
all_dat_df = all_df_query(search_string, save_csv= True)
print(all_dat_df.shape)

# Testing out CLUSTALOMEGA functions
##################################################################

# Post a job to Clustal Omega
job_id = post_alignment_request(fasta_file)
print("Clustal Omega job successfully submitted with Job ID: " + job_id)

# Wait for job to finish running, then print it
clustal_omga_file = clustalo_alignment(job_id) 
while clustal_omga_file == False:
    time.sleep(3)
    clustal_omga_file = clustalo_alignment(job_id) 

print(clustal_omga_file)
