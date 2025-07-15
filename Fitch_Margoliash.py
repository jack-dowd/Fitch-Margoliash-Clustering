#!/usr/bin/env python3
"""
================================================================================
Title          : Fitch_Margoliash.py
Description    : Multiple Sequence Alignment
Author         : Jack Dowd (R11690912)
Date           : November 7, 2024
Version        : 1.0
Usage          : python3 Fitch_Margoliash.py -i input.fna -o output.fna -s BLOSUM62.mtx
Notes          : This tool requires a substitution matrix file (Ex: BLOSUM62).
Python Version : 3.10.12
================================================================================
"""

import argparse
import math
import os
import copy

class Sequence:
    def __init__(self, identifier, sequence_score, sequence_data, sequence_length, index):
        self.identifier = identifier
        self.sequence_score = sequence_score
        self.sequence_data = sequence_data.replace("\n", "").replace("\r", "")
        self.sequence_length = sequence_length
        self.index = index
    
class SeqPair:
    def __init__(self, seq1, seq2, similarity_score, distance):
        self.seq1 = seq1
        self.seq2 = seq2
        self.similarity_score = similarity_score
        self.distance = distance
        
# get Sequence
def getSequence(identifier, sequence_data, index):
    sequence_length = len(sequence_data) - (str.count(sequence_data, '\n') + str.count(sequence_data, '\r'))
    return Sequence(identifier, 0, sequence_data, sequence_length, index)
def getSequenceFromInputList(input_list, index):
    string = input_list[index]
    split = str.split(string, '\n', 1)
    # identifier
    identifier = split[0].replace("\r", "")
    # sequence data
    if(len(split) >= 2):
        sequence_data = split[1]
    else:
        sequence_data = ""
    # sequence length
    sequence_length = len(sequence_data) - (str.count(sequence_data, '\n') + str.count(sequence_data, '\r'))
    # return Sequence
    return Sequence(identifier, 0, sequence_data, sequence_length, index)
        
# get SeqPair
def getSeqPair(seq1, seq2, score):
    return SeqPair(seq1, seq2, score, 0)

# print a sequence
def printSequence(sequence):
    string = ">" + sequence.identifier + "; score=" + str(sequence.sequence_score) + "\n" + sequence.sequence_data + "\n"
    print(string)

# to check if the output directory is the current directory or not
def isCurrentDirectory(output_file):
    numSlashes = str.count(output_file, "/")
    # if the format is "/output.fna" (don't count the slash)
    if(numSlashes == 1 and output_file[0] == "/"):
        return True
    # if the format is "/path/to/output.fna" or "path/to/output.fna"
    elif(numSlashes > 0):
        return False
    # if the format is "output.fna"
    elif (numSlashes == 0):
        return True
    
# check if the output directory exists
def outputDirExists(output_file):
    # if the output directory is the current directory
    if(isCurrentDirectory(output_file) == True):
        return True
    # get directory name
    dir_name = os.path.dirname(output_file)
    # check if directory name exists
    if(os.path.exists(dir_name)):
        return True
    else:
        return False

# check if the input file exists
def inputFileExists(input_file):
    if(os.path.isfile(input_file)):
        return True
    else:
        return False

# return the number of sequences in the input file
def getNumSeqs(input_file):
    num_seqs = 0
    with open(input_file, 'r', newline='\n') as file:
        for seq in file:
            if(seq.startswith('>')):
                num_seqs += 1
    return num_seqs

# create output file
def generateOutputFile(output_file, output_string):
    with open(output_file, 'w', newline='\n') as file:
        file.write(output_string)

# sort sequence list from highest length to lowest length
def sortSequenceList(sequence_list):
    return sorted(sequence_list, key=lambda sequence: sequence.sequence_length)
    
# parse score matrix and generate lookup dict
def parse_score_matrix(score_matrix_string):
    rows = score_matrix_string.strip().splitlines() # split the text into lines and remove leading/trailing whitespace
    headers = rows[0].split()
    
    lookup_dict = {}
    for row in rows[1:]: # skip the header row
        values = row.split()
        row_key = values[0] # row identifier (A, R, N, D,...)
        lookup_dict[row_key] = {headers[i]: int(values[i+1]) for i in range(len(headers))}
    
    # add 'X' if it doesn't already exist
    if list.count(headers, 'X') == 0:
        # add 'X' to dict (always -1)
        lookup_dict['X'] = {header: -1 for header in headers}
        lookup_dict['X']['X'] = -1
        
        # update each existing row to include 'X' as -1
        for key in lookup_dict:
            lookup_dict[key]['X'] = -1

        # add gap penalties at ['X']['-']
        lookup_dict['X']['-'] = lookup_dict['-'][headers[0]]
        lookup_dict['-']['X'] = lookup_dict['-'][headers[0]]
        
        # add 'X' to headers
        headers.append('X')

    return lookup_dict

# parse score matrix to a 1d array (doesn't include '-' scores)
def parse_score_matrix_to_1d_array(score_matrix_string):
    rows = score_matrix_string.strip().splitlines() # split the text into lines and remove leading/trailing whitespace
    score_array = []

    for row in rows[1:-1]: # exclude header and gap penalties (because S_max and S_min exclude gap penalties)
        scores = list(map(int, row.split()[1:-1]))
        score_array.extend(scores)

    return score_array

# find the score value from the score matrix
def lookup(row, col, lookup_dict):
    if row in lookup_dict and col in lookup_dict[row]:
        return lookup_dict[row][col]
    else:
        return None

def calculateScore(aligned_sequences, lookup_dict):
    score = 0
    aligned_strings = []
    for sequence in aligned_sequences:
        aligned_strings.append(sequence.sequence_data)
        
    for col in range(len(aligned_strings[0])):
        column_chars = [seq[col] for seq in aligned_strings]
        for i, char1 in enumerate(column_chars):
            for char2 in column_chars[i+1:]:
                # Penalize only non-shared gaps
                score += lookup(char1, char2, lookup_dict)
                
    return score
        
# align sequences
def needleman_wunsch(seq1, seq2, lookup_dict, replace_gaps_with_X):

    # change gaps to '-' (used during Feng-Doolittle)
    gapChar = '-'
    if replace_gaps_with_X == True:
        gapChar = 'X'
    
    m = len(seq1)
    n = len(seq2)

    # create score matrix
    score_matrix = [[0 for j in range(n+1)] for i in range(m+1)]

    # initialize the first row and column with gap penalties from lookup_dict
    gap_penalty = lookup_dict['-']['X']
    for i in range(1, m+1):
        score_matrix[i][0] = score_matrix[i-1][0] + gap_penalty
    for j in range(1, n+1):
        score_matrix[0][j] = score_matrix[0][j-1] + gap_penalty

    # fill the score matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_mismatch = score_matrix[i-1][j-1] + lookup_dict[seq1[i-1]][seq2[j-1]]
            insertion_1    = score_matrix[i-1][j]   + gap_penalty
            insertion_2    = score_matrix[i][j-1]   + gap_penalty
            score_matrix[i][j] = max(match_mismatch, insertion_1, insertion_2)

    # traceback to get the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = m
    j = n
    
    while(i > 0 and j > 0):
        current_score = score_matrix[i][j]
        if current_score == score_matrix[i-1][j-1] + lookup_dict[seq1[i-1]][seq2[j-1]]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif current_score == score_matrix[i-1][j] + gap_penalty:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = gapChar + aligned_seq2
            i -= 1
        elif current_score == score_matrix[i][j-1] + gap_penalty:
            aligned_seq1 = gapChar + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
    # Add padding if stopping early leaves lengths unequal
    while i > 0:
        aligned_seq1 = seq1[i-1] + aligned_seq1
        aligned_seq2 = gapChar + aligned_seq2
        i -= 1
    while j > 0:
        aligned_seq1 = gapChar + aligned_seq1
        aligned_seq2 = seq2[j-1] + aligned_seq2
        j -= 1
    
    # return aligned sequences
    return aligned_seq1, aligned_seq2, score_matrix[m][n]

# all vs all global alignment
def all_vs_all_global_alignment(sequence_list, lookup_dict):
    
    # num of sequences
    n = len(sequence_list)
    
    # iterate over all pairs of sequences
    score = 0
    aligned_pairs = []
    for i in range(n):
        for j in range(i+1, n):
            sequences = copy.deepcopy(sequence_list)
            seq1 = getSequence(sequences[i].identifier, sequences[i].sequence_data, sequences[i].index)
            seq2 = getSequence(sequences[j].identifier, sequences[j].sequence_data, sequences[j].index)
            seq1.sequence_data, seq2.sequence_data, s = needleman_wunsch(seq1.sequence_data, seq2.sequence_data, lookup_dict, replace_gaps_with_X=False)
            seq1.sequence_score, seq2.sequence_score = s, s
            aligned_pairs.append(getSeqPair(seq1, seq2, s))
            score += s
            
    return aligned_pairs, score

# distance matrix for pairs
def compute_distance_matrix(seq_pairs, score_matrix_1d_array, num_seqs):
    
    min_value = float(min(score_matrix_1d_array)) # the smallest value in the score matrix
    max_value = float(max(score_matrix_1d_array)) # the largest value in the score matrix

    # set k
    k = 10.0

    # compute distances and normalize similarity scores
    for i in range(len(seq_pairs)):
        length = float(seq_pairs[i].seq1.sequence_length)
        score = float(seq_pairs[i].similarity_score)
        S_min = float(length * min_value)
        S_max = float(length * max_value)
        S_n = float((score - S_min) / (S_max - S_min))
        seq_pairs[i].distance = int(-k * math.log(S_n))

    distance_matrix = [[0.0 for j in range(num_seqs)] for i in range(num_seqs)]
    for i in range(len(seq_pairs)):
        distance_matrix[seq_pairs[i].seq1.index][seq_pairs[i].seq2.index] = seq_pairs[i].distance
        
    return seq_pairs, distance_matrix


# Shorten a nested binary list to a string
# Ex: "[0, [1, [2, 3]]]" -> "[0[1[23]]]"
def joinpair(pair):
    if isinstance(pair, int):
        return pair
    if is_1d_list(pair):
        return "[" + str(pair[0]) + str(pair[1]) + "]"
    else:
        l = pair[0]
        r = pair[1]
        l = joinpair(l)
        r = joinpair(r)
        return "[" + str(l) + str(r) + "]"



def print_distance_matrix(distance_matrix, labels):
    print()
    lbl = []
    for i in range(len(labels)):
        if isinstance(labels[i], list):
            lbl.append("[" + str(labels[i][0]) + str(labels[i][1]) + "]")
        else:
            lbl.append("[" + str(labels[i]) + "]")
            
    for i in range(len(distance_matrix)):
        print(f"\t{lbl[i]}",end='')
    
    print()
    for i in range(len(distance_matrix)):
        print(f"{lbl[i]}",end='')
        for j in range(len(distance_matrix[i])):
            print(f"\t{distance_matrix[i][j]}",end='')
        print()
        
def print_distance_dict(distance_dict):
    for i in distance_dict:
        print(i, distance_dict[str(i)])




def mirror_matrix(matrix):
    # Get the size of the matrix
    n = len(matrix)
    
    # Iterate over the upper triangle of the matrix (excluding the diagonal)
    for i in range(n):
        for j in range(i + 1, n):
            # Swap elements (i, j) and (j, i)
            matrix[j][i] = matrix[i][j]
    
    return matrix


def list_to_2d_dict(labels, integer_list):
    # Create a dictionary with each label mapped to another dictionary
    return {
        str(labels1): {str(labels[index]): value for index, value in enumerate(labels2)} for labels1, labels2 in zip(labels, integer_list)
    }
    
def print_2d_dict(dic):
    print(dic)

# Fitch-Margoliash
def fitch_margoliash_clustering(distance_matrix, labels):
    
    n = len(distance_matrix)
    final_labels = []
    a = 0
    while n > 1:
        # Step 1: Find the closest sequences
        min_dist = math.inf
        x, y = -1, -1
        
        
        distance_dict = list_to_2d_dict(labels, distance_matrix)
        #print_2d_dict(distance_dict)
        
        for i in range(n):
            for j in range(i + 1, n):
                if distance_dict[str(labels[i])][str(labels[j])] < min_dist:
                    min_dist = distance_dict[str(labels[i])][str(labels[j])]
                    x, y = labels[i], labels[j]



        # Step 2: Create a new label for the merged sequences
        new_label = [x, y]
        final_labels.append(new_label)
        
        # Update labels by removing old and appending new
        new_labels = [labels[i] for i in range(n) if labels[i] != x and labels[i] != y]
        new_labels.append(new_label)
        #if isinstance(new_labels[0], int):
        #print(new_labels, "|", final_labels)
        #print("-----------------------------")
        #print(new_labels)

        # Step 3: Generate a new matrix merging x & y
        new_distance_matrix = [[0] * (n - 1) for _ in range(n - 1)]
        #print()
        #print()
        #print(new_labels_ints)
        for i in range(n-1):
            for j in range(i+1, n-2):  
                #print(new_labels_ints[i], new_labels_ints[j])
                new_distance_matrix[i][j] = distance_dict[str(new_labels[i])][str(new_labels[j])]
        
        for i in range(len(new_distance_matrix)-1):
        #    print(">>>>", i)
            total = 0
        #    print(i, x, y)
            total = (distance_dict[str(x)][str(new_labels[i])] + distance_dict[str(y)][str(new_labels[i])]) / 2
            new_distance_matrix[i][n-2] = total

        #print_distance_matrix(new_distance_matrix, new_labels)
                
        #print()
        lbls = []
        for i in range(len(new_labels)):
            if isinstance(new_labels[i], list):
                lbls.append(joinpair(new_labels[i]))
            else:
                lbls.append("[" + str(new_labels[i]) + "]")
                
        #for i in range(len(new_distance_matrix)):
            #print(f"\t{lbls[i]}",end='')
        
        #print()
        #for i in range(len(new_distance_matrix)):
        #    if len(lbls[i]) > 7:
        #        print(f"{lbls[i][0:7]}\t",end='')
        #    else:
        #        print(f"{lbls[i]}\t",end='')
        #    for j in range(len(new_distance_matrix[i])):
        #        print(f"{new_distance_matrix[i][j]}\t",end='')
        #    print()
            
        # Update distance matrix and labels for next iteration
        distance_matrix = new_distance_matrix
        labels = new_labels
        n -= 1
        a += 1
    #print()
    #print(final_labels)
    return final_labels


# get the index of a label within a list of labels
def get_label_index(labels, label):
    for i in range(len(labels)):
        if(labels[i] == label):
            return i
        else:
            return None
    return None

# check if the argument is a 1D list
def is_1d_list(arg):
    if isinstance(arg, list):  # check if arg is a list
        return all(not isinstance(item, list) for item in arg)  # ensure no item is a list
    return False  # not a list

def is_integer(value):
    if isinstance(value, int):
        return True
    elif isinstance(value, float):
        return value.is_integer()
    return False

def flatten(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten(item))  # Recursively flatten the item
        else:
            flat_list.append(item)  # Add the item directly
    return flat_list

# perform MSA from guide_tree
def msa(sequences, node, lookup_dict):

    if is_integer(node):
        node = [node]
        
    # If `node` is in the form [0,1]
    if is_1d_list(node) and len(node) == 2:
        
        # get sequences from guide tree values
        seq1_index = node[0]
        seq2_index = node[1]
        seq1 = sequences[node[0]]
        seq2 = sequences[node[1]]
        
        aligned_seq1, aligned_seq2, _ = needleman_wunsch(seq1, seq2, lookup_dict, replace_gaps_with_X=False)
        
        # Create a list to hold the aligned sequences in their correct indices
        aligned_sequences = [0 for i in range(len(sequences))]
        aligned_sequences[seq1_index] = aligned_seq1
        aligned_sequences[seq2_index] = aligned_seq2
        
        return aligned_sequences, [seq1_index, seq2_index]
    
    # If the node is a list with a length of 1
    elif len(node) == 1:

        # fill in the aligned sequences for the left subtree (node[0])
        final_aligned_sequences = [0 for i in range(len(sequences))]
        final_aligned_sequences[node[0]] = sequences[node[0]]

        return final_aligned_sequences, node
    
    # If the node is a list but not a 1-dimensional list, like [0, [1, 2]] or [[0, 1], [2, 3]]
    else:
        
        # get subtrees
        left_subtree = node[0]
        right_subtree = node[1]
        
        # align the left and right subtrees recursively
        left_aligned, left_subtree = msa(sequences, left_subtree, lookup_dict)
        right_aligned, right_subtree = msa(sequences, right_subtree, lookup_dict)
        
        # Create a list to hold the aligned sequences
        final_aligned_sequences = [0 for i in range(len(sequences))]
        
        # Case: "[[0, 1], [2, 3]]"
        if len(left_subtree) == 2 and len(right_subtree) == 2:
            # We will be comparing each index in the left subtree to *all* indexes in the original right subtree (node[1])
            
            # Flatten the right subtree
            right_subtree_flattened = flatten(node[1])

            # First half of the case
            seq1 = ''
            seq2 = ''
            max_score = -math.inf
            max_score_index = -1
            left_index = -1
            
            # compare each sequence in left_subtree with the sequences in right_subtree
            for i in range(len(left_subtree)): # this length is always 2
                for j in range(len(right_subtree_flattened)):
                    a, b, s = needleman_wunsch(left_aligned[left_subtree[i]], right_aligned[right_subtree_flattened[j]], lookup_dict, replace_gaps_with_X=True)
                    if s > max_score:
                        seq1 = a
                        seq2 = b
                        max_score = s
                        max_score_index = right_subtree_flattened[j]
                        left_index = left_subtree[i]

            # Get right subtree sequences to align
            leftover_idxs = []
            leftover_seqs = []
            for i in range(len(right_subtree_flattened)):
                if max_score_index != right_subtree_flattened[i]:
                    leftover_idxs.append(right_subtree_flattened[i])
            for i in range(len(leftover_idxs)):
                leftover_seqs.append(right_aligned[leftover_idxs[i]])
            leftover_left_idx = left_subtree[0]
            
            # Get left subtree sequence to align
            if left_index == left_subtree[0]:
                leftover_left_idx = left_subtree[1]
            else:
                leftover_left_idx = left_subtree[0]
            leftover_left_seq = left_aligned[leftover_left_idx]
            
            # Align gaps in left subtree
            left_aligned[left_index] = seq1.replace('X', '-')
            left_aligned[leftover_left_idx] = align_gaps(seq1, leftover_left_seq).replace('X', '-')
            
            # Align gaps in right subtree
            for i in range(len(leftover_seqs)):
                right_aligned[leftover_idxs[i]] = align_gaps(seq2, leftover_seqs[i]).replace('X', '-')
            right_aligned[max_score_index] = seq2.replace('X', '-')
            
            # Update left and rigth subtrees maximum score indexes
            right_subtree = [max_score_index]
            left_subtree = [left_index]
            
        # Case: "[[0], [1, 2]]"
        elif len(left_subtree) == 1 and len(right_subtree) == 2:
            left_subtree_flattened = [node[0]]
            right_subtree_flattened = flatten(node[1])
            
            seq1 = ''
            seq2 = ''
            max_score = -math.inf
            max_score_index = -1
            for index in right_subtree_flattened:
                a, b, s = needleman_wunsch(left_aligned[left_subtree[0]], right_aligned[index], lookup_dict, replace_gaps_with_X=True)
                if s > max_score:
                    seq1 = a
                    seq2 = b
                    max_score = s
                    max_score_index = index
            left_index = left_subtree[0]

            # Get right subtree sequences to align
            leftover_idxs = []
            leftover_seqs = []
            for i in range(len(right_subtree_flattened)):
                if max_score_index != right_subtree_flattened[i]:
                    leftover_idxs.append(right_subtree_flattened[i])
            
            # Align gaps in right subtree
            for i in range(len(leftover_idxs)):
                leftover_seqs.append(right_aligned[leftover_idxs[i]])
            for i in range(len(leftover_seqs)):
                right_aligned[leftover_idxs[i]] = align_gaps(seq2, leftover_seqs[i]).replace('X', '-')
            left_aligned[left_index] = seq1.replace('X', '-')
            right_aligned[max_score_index] = seq2.replace('X', '-')
            
            # Update right_subtree to max_score_index
            right_subtree = [max_score_index]
        
        # fill in final_aligned_sequences
        for i in range(len(final_aligned_sequences)):
            if left_aligned[i] != 0 and right_aligned[i] != 0:
                raise TypeError("left_aligned[i] and right_aligned[i] should not be strings. At least one of them should be the value 0.")
            elif left_aligned[i] != 0:
                final_aligned_sequences[i] = left_aligned[i]
            elif right_aligned[i] != 0:
                final_aligned_sequences[i] = right_aligned[i]

        # return the final alignment and the updated subtrees, which should be in the form [0,1] no matter what
        return final_aligned_sequences, [left_subtree[0], right_subtree[0]]

# insert a character into a string at a specific index
def insert_char_at_index(original_string, char_to_insert, index):
    if index < 0:
        index = 0
    elif index > len(original_string):
        index = len(original_string)
    new_string = original_string[:index] + char_to_insert + original_string[index:]
    return new_string

# align gaps for Feng-Doolittle
def align_gaps(seq_with_gaps, seq_to_align):
    aligned_seq = ""
    gaps_added = 0
    for i in range(len(seq_with_gaps)):
        if seq_with_gaps[i] == 'X':
            aligned_seq += '-'
            gaps_added += 1
        else:
            aligned_seq += seq_to_align[i - gaps_added]
    return aligned_seq


def pop_char(string, index):
    return string[:index] + string[index + 1:], string[index]



def max_length(sequences):
    max_len = 0
    for sequence in sequences:
        if sequence != 0:
            if len(sequence) > max_len:
                max_len = len(sequence)
    return max_len

def min_length(sequences):
    min_len = math.inf
    for sequence in sequences:
        if len(sequence) < min_len:
            min_len = len(sequence)
    return min_len
        













def main():   
    # print R number
    print("Assignment 4 :: R11690912")
    
    # arg parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-s', '--score_matrix_file', type=str, required=True)
    parser.add_argument('-n', '--number_of_iterations', type=int, required=False)
    
    # get args
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    score_matrix_file = args.score_matrix_file
    number_of_iterations = args.number_of_iterations

    # check if input file exists
    if inputFileExists(input_file) == False:
        raise ValueError("Input file not found.")
    # check if output file directory exists
    if outputDirExists(output_file) == False:
        raise ValueError("Output directory not found.")
    # check if score matrix file exists
    if inputFileExists(score_matrix_file) == False:
        raise ValueError("Score matrix file not found.")
    # if number of iterations is None, set to a default value (1)
    if number_of_iterations == None:
        number_of_iterations = 1
    
    
    
    
    
    
    
    # get number of sequences in input file
    num_seqs = getNumSeqs(input_file)
    
    # create a string of the score matrix file text
    score_matrix_string = ""
    with open(score_matrix_file, 'r', newline='\n') as file:
        while True:
            line = file.readline()
            if not line:
                break # EOF reached
            if str.isspace(line) == False: # (ignore blank lines)
                score_matrix_string += line
    
    # parse the score matrix and get lookup dict
    lookup_dict = parse_score_matrix(score_matrix_string)   
    
    # create a string of the input file text
    input_string = ""
    with open(input_file, 'r', newline='\n') as file:
        file.read(1) # skip the first '>' (to avoid creating a blank Sequence instance later at sequence_list[0])
        while True:
            line = file.readline()
            if(not line):
                break # EOF reached
            input_string += line
    
    # split entire string into a list of Sequence instances (identifier, description, sequence data, sequence length)
    input_list = str.split(input_string, ">")
    sequence_list = [0 for i in range(num_seqs)]
    for i in range(num_seqs):
        sequence_list[i] = getSequenceFromInputList(input_list, i)
        
    # begin iterations loop
    prev_msa = [0 for i in range(num_seqs)]
    for n in range(max(number_of_iterations,1)):
            
        # load previous iteration's MSA to current iteration's starting sequences
        if n != 0:
            for i in range(num_seqs):
                sequence_list[i].sequence_data = prev_msa[i]
        
        # all vs all alignment
        seqs = [0 for i in range(num_seqs)]
        for i in range(num_seqs):
            seqs[i] = sequence_list[i]
        seq_pairs, pairs_score = all_vs_all_global_alignment(seqs, lookup_dict) # seq_pairs is a list of SeqPair instances
        
        
        # get score matrix values as a 1d array 
        #   (used to find S_min)
        score_matrix_1d_array = parse_score_matrix_to_1d_array(score_matrix_string)
        
        # compute distance matrix (distances are also saved in `seq_pairs[i].distance`)
        seq_pairs, distance_matrix = compute_distance_matrix(seq_pairs, score_matrix_1d_array, num_seqs)



        # Fitch-Margoliash
        labels = [i for i in range(num_seqs)]
        guide_tree = fitch_margoliash_clustering(distance_matrix, labels) # returns a list in the format [[0,1],[2,3]] or similar at the last index (guide_tree[-1])
            
        string_list = [0 for i in range(num_seqs)]    
        for i in range(num_seqs):
            string_list[i] = sequence_list[i].sequence_data
        
        # Feng-Doolittle / MSA
        msa_strings, subtrees = msa(string_list, guide_tree[-1], lookup_dict)


        # Add padding just in case
        max_len = max_length(msa_strings)
        for i in range(num_seqs):
            if len(msa_strings[i]) < max_len:
                for j in range(max_len - len(msa_strings[i])):
                    msa_strings[i] = msa_strings[i] + '-'


        # Remove columns with all gaps
        for i in range(num_seqs):
            msa_strings[i] = str.replace(msa_strings[i], 'X', '-')
        all_gaps_indexes = []
        for index in range(len(msa_strings[0])):
            all_gaps = True
            for seq in msa_strings:
                if seq[index] != '-':
                    all_gaps = False
            if all_gaps:
                all_gaps_indexes.append(index)
        
        for i in range(num_seqs):
            for index in reversed(all_gaps_indexes):
                msa_strings[i], popped_char = pop_char(msa_strings[i], index)
                    
        # calculate score
        score = 0
        for col in range(len(msa_strings[0])):
            s = 0
            for i in range(num_seqs):
                for j in range(i+1, num_seqs):
                    s += lookup(msa_strings[i][col], msa_strings[j][col], lookup_dict)
            score += s  

        # create output string
        output_string = ""
        for i in range(num_seqs):
            sequence = sequence_list[i]
            sequence_string = ">" + sequence.identifier + "; score=" + str(score) + "\n" + msa_strings[i] + "\n"
            output_string += sequence_string
        
        # print output string
        print(output_string)
        
        # update sequences for next iteration
        for i in range(len(msa_strings)):
            prev_msa[i] = msa_strings[i]
            

    # generate output file
    generateOutputFile(output_file, output_string)

# value error handling
if __name__ == '__main__':
    try:
        main()
    except ValueError as e:
        print("ValueError:", e)
