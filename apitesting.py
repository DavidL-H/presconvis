from uniprotAPI import all_fasta_query, uniprot_query, all_df_query, response_to_df
import pandas as pd

# Test out fetching of uniprot protein sequences from custom function:
search_string = "uniref_cluster_90:UniRef90_Q8X825"

# Print query metadata
response = uniprot_query(search_string)
print(response.headers)

# Write query to fasta file
fasta_file = all_fasta_query(search_string, force=False)
print(fasta_file)

# Extract key data to dataframe:


df_temp = response_to_df(uniprot_query(search_string))
page2 = "uniref_cluster_90:UniRef90_Q8X825&cursor=1mkycb2xwxboutolnex194jt6msouvkwiblg&size=25"
df_temp2 = response_to_df(uniprot_query(page2))

frames = [df_temp,df_temp2]
query_df = pd.concat(frames)
#all_df = all_df_query(search_string)


# Write json for manual inspection
#with open("sample.json", "w") as outfile:
    #json.dump(json_file, outfile)




