'''
Perform analysis on a Multiple sequence alignment in the form of a .clustal_num file
David Lennox-Hvenekilde
220208
'''
import numpy as np
import pandas as pd
root = "./endpoints/"

example = "clustalo-R20220801-162757-0097-51791681-p1m.clustal_num"
with open(root+example, "r") as alignment:
    first_entry = "\n"
    while first_entry == '\n' or first_entry == 'CLUSTAL':
        first_line = alignment.readline()

        all_split = first_line.split(" ")
        print(all_split[len(all_split)-1])

        first_entry = first_line.split(" ")[0]

print(first_entry)

# Pseudocode: ######################################################

# Generate a main_matrix with n row (items in MSA) and no columns
# Need to loop through the document
    # Identify lines starting with the uniprot id: "base_id"
        # When based_id is indentified, generate an empty sub_matrix
        # loop through the next lines until a space:
        # Extract the content between spaces and append row to the sub_matrix
    # Loop through the columns of the matrix, 
        # if element of first row != "-"
            #  append column to main_matrix


# Check number of rows in file
with open(root+example, "r") as alignment:
    line_count = len(alignment.readlines())

# Create empty array
A = np.array([None] * 61)
trigger = True

# Open file
with open(root+example, "r") as alignment:
    first_entry = "\n"
    # Read each row
    for i in range(line_count):
        line = alignment.readline()
        line_split =  line.split(" ")
        
        # Qualifiers for containing relevant information
        if line_split[0] != '\n' and line_split[0] != 'CLUSTAL' and line_split[0] != '':
            id_name = line_split[0]
            # Grab the body of the alingment
            body = line_split[len(line_split)-1]
            # Split again by tab
            second_split = body.split("\t")
            amino_acids = second_split[0]

            # Concatenate name and aminoacids, and append to matrix as row
            # This breaks for the last part of the MSA, that is not the same length as the other
            # For this last part, we need to generate a new matrix of a difference column length
            try:
                A = np.vstack([A, [id_name]+list(amino_acids)])
            except:
                if trigger:
                    A2 = np.array([None] * len([id_name]+list(amino_acids)))
                    trigger = False
                else:
                    A2 = np.vstack([A2, [id_name]+list(amino_acids)])

# Can use numpy hstack to concatenate matrices horizontally
A.shape
A2.shape

# Count unique sequences
seq_len = len(set(A[:,0]))

# Loop through the matrix and reshape it based on the unqiue sequence identifiers
new_a = np.hstack((A,A))

