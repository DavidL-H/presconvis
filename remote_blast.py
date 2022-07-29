# Functions for Running BLAST jobs remotely with Bio.Blast
# https://www.tutorialspoint.com/biopython/biopython_overview_of_blast.htm

from Bio.Blast import NCBIWWW
import time

with open("UniRef90_Q8X825.fasta", "r") as fasta_file:
    fasta_file.readlines()

seq_test = "MAHRPRWTLSQVTELFEKPLLDLLFEAQQVHRQHFDPRQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLESARKAKAAGSTRFCMGAAWKNPHERDMPYLEQMVQGVKAMGLEACMTLGTLSESQAQRLANAGLDYYNHNLDTSPEFYGYIITTRTYQERLDTLEKVREAGIKVCSGGIVGLGETVKDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPEEDKDLQLFRKLGLNPQQTAVLAGDNEQQQRLEQALMTPDTDEYYNAAAL"

# This takes ages....
start = time.time()
result_handle = NCBIWWW.qblast("blastp", "nt", seq_test) 
end = time.time()
print("BLAST query took: "+ str(end-start))

# Save results to xml
with open('./endpoints/results.xml', 'w') as save_file:
    blast_results = result_handle.read() 
    save_file.write(blast_results)