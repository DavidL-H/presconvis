from uniprotAPI import all_fasta_query, uniprot_query

search_string = "uniref_cluster_90:UniRef90_Q8X825"
fasta_file = all_fasta_query(search_string, force=True)
print(fasta_file)

