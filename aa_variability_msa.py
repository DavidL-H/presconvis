'''
Perform analysis on a Multiple sequence alignment in the form of a .clustal_num file
David Lennox-Hvenekilde
220208
'''

from math import log
from typing import final
import numpy as np
import pandas as pd
import os
from collections import Counter
import blosum as bl

root = "./endpoints/"

# Turn the clustal omega alignment into a matrix/dataframe
def clustalo_to_matrix(clustalo_file, save_matrix = True, force = False):
    '''
    This funciton take a Clustal Omega multiple sequence alignment file (.clustal_num)
    and converts it a large numpy matrix, with one row per sequence id and one column per sequence position.
    Saves the matrix as a .csv file, if wanted.
    '''
    alignment_name = clustalo_file.split(".")[0]
    # Check existing file
    if os.path.isfile(root + alignment_name + "_MSA" +".csv"):
        if force:
            os.remove(root + alignment_name + "_MSA" +".csv")
            print("Existing clustal omega file deleted, and remade")
        elif not force:
            print("Clustal omega matrix file already exists")
            return pd.read_csv(root + alignment_name+ "_MSA" +".csv")

    # Check number of rows in file
    with open(root+clustalo_file, "r") as alignment:
        line_count = len(alignment.readlines())

    # Create empty array
    A = np.array([None] * 61)
    trigger = True

    # Open file
    with open(root+clustalo_file, "r") as alignment:
        # Read each row
        for _ in range(line_count):
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
    # Delete first row, used for initiating the array
    A = np.delete(A, (0), axis=0)

     # Count unique sequences
    seq_len = len(set(A[:,0]))

    # Number of times
    start = 0 
    end = seq_len
    new_array = A[start:end,:]
    for _ in range(int(len(A)/seq_len)-1):
        start += seq_len
        end += seq_len
        temp_array = A[(start):end,1:]

        # concatenate matrices horizontally
        new_array = np.hstack((new_array,temp_array))

    # Saving the dataframe
    if save_matrix:
        print("Matrix saved as: " + root + alignment_name + "_MSA" +".csv")
        pd.DataFrame(new_array).to_csv(root + alignment_name + "_MSA" +".csv")
        return pd.DataFrame(new_array)
    else:
        return pd.DataFrame(new_array)


# Reformat the matrix/dataframe
def clustalo_df_formating(clustalo_df, query_seq = "default"):
    '''
    Takes the dataframe/matrix output from clustalo_to_matrix() and reformats it
    to root the sequence in the "query_seq". By "default" this is first sequence
    to show up in the original clustal omega alignment, but it can be specified as well
    '''     
    clustalo_df = clustalo_df.set_index("0")
    del clustalo_df["Unnamed: 0"]

    # choose row to normalize/root by
    if query_seq == "default":
        query_seq = clustalo_df.index[0]

    temp_df = clustalo_df.loc[[query_seq]]

    # Drop empty columns, in regards to query_seq
    column_drop = []
    for column in range(len(temp_df.columns.values)):
        if temp_df.iloc[0][column] == "-":
            column_drop.append(column)

    clustalo_df_normalized = clustalo_df.drop(clustalo_df.columns[column_drop],axis = 1)
    # Recount columns in regard to query seq
    clustalo_df_normalized.columns = list(range(1,len(clustalo_df_normalized.columns)+1))

    # We now have our properly formatted dataframe
    return clustalo_df_normalized


# Generate the final formatted data frame
def final_MSA_df(formatted_df, file_name = "Final_Score_df", save_df = True, force = False):
    '''
    Format the data frame to give a better overview of amino acid occurance at each position
    and calculate a variability score at each said position
    '''

    # Check existing file
    if os.path.isfile(root + file_name +".csv"):
        if force:
            os.remove(root + file_name +".csv")
            print("Existing data frame file deleted, and remade")
        elif not force:
            print("Data frame file already exists. Return saved dataframe")
            return pd.read_csv(root + file_name +".csv")


    # Lets do some summary statistics on the dataframe
    # create a new data frame, initialized with AA index of query sequence
    final_df = pd.DataFrame(formatted_df.iloc[0].values.tolist(), columns=["Residue"])

    # Calculate residue variability scores
    blosum_matrix = bl.BLOSUM(62)

    # Count residues for each position
    all_scores = []
    aas = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-"]
    list_of_postiion_values = []

    for pos in formatted_df.columns.tolist():
        # Count occurance of residue at position: pos
        pos_dict = dict(Counter(formatted_df[pos].tolist()))

        # Split occurances into a nested list for generating a df
        aas_n = [0] * len(aas)
        dict_keys = list(pos_dict.keys()) 
        for k_i in dict_keys:
            aas.index(k_i)
            aas_n[aas.index(k_i)] = pos_dict[k_i]

        list_of_postiion_values.append(aas_n)

        # Delete empty alignment positions, as this doesn't play well with BLOSUM 62 matrix
        del pos_dict["-"]

        # Generate a score based on BLOSUM 62 and prevalence
        total_score = 0
        for key in pos_dict:
            for key_2 in pos_dict:
                total_score += blosum_matrix[key+key_2]*(pos_dict[key]+pos_dict[key_2])

        total_score = round(total_score/len(formatted_df),2)

        all_scores.append(total_score)

    # Add counts of residue occurance at each residue position, where each residue option is a column
    # and each row corresponds to the residue position
    df_counts_positions = pd.DataFrame(list_of_postiion_values, columns = aas)

    # Normalize score by number of proteins in alignment
    final_df["Scores"] = all_scores
    final_df["Position"] = final_df.index.tolist()

    # Merge scores and counts data frames
    final_df = pd.concat([final_df, df_counts_positions], axis=1)
    # Save for inspection
    if save_df:
        final_df.to_csv(root + file_name + ".csv")
        return pd.read_csv(root + file_name +".csv")
    else:
        return final_df

