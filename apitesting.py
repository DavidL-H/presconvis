from uniprotAPI import all_fasta_query, uniprot_query, all_df_query
import clustalomegaAPI as clust
import pandas as pd

# Testing out UNIPROT functions
##################################################################

# Test out fetching of uniprot protein sequences from custom function:
search_string = "uniref_cluster_90:UniRef90_Q8X825"

# Print query metadata
response = uniprot_query(search_string)
print(response.headers)

# Write query to fasta file
fasta_file = all_fasta_query(search_string, force=False)
print(fasta_file)

# Extract key data to dataframe:
search_string = "uniref_cluster_90:UniRef90_Q8X825"
all_dat_df = all_df_query(search_string)
print(all_dat_df.shape)

# Testing out CLUSTALOMEGA functions
##################################################################


