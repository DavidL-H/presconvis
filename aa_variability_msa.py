'''
Perform analysis on a Multiple sequence alignment in the form of a .clustal_num file
David Lennox-Hvenekilde
220208
'''
import pandas as pd
root = "./endpoints/"

example = "clustalo-R20220801-162757-0097-51791681-p1m.clustal_num"
with open(root+example, "r") as alignment:
    for _ in range(4):
        first_line = alignment.readline()

base_id = first_line.split(" ")[0]

# Pseudocode: #

# Generate a main_matrix with n row (items in MSA) and no columns
# Need to loop through the document
    # Identify lines starting with the uniprot id: "base_id"
        # When based_id is indentified, generate an empty sub_matrix
        # loop through the next lines until a space:
        # Extract the content between spaces and append row to the sub_matrix
    # Loop through the columns of the matrix, 
        # if element of first row != "-"
            #  append column to main_matrix


with open(root+example, "r") as alignment:
    i = 0
    while id != base_id:
        i = i + 1
        line = alignment.readline()
        id = line.split(" ")[0]
    print("Next line is: "+str(i)+" ",line)