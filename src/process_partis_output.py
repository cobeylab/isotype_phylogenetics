#!/usr/bin/python

"""Processes partis output to produce separate fasta files with the sequences of each clone, along with inferred naive sequence for the clone
"""
import sys
import csv
from heapq import nlargest
import re
csv.field_size_limit(sys.maxsize)

# Process the n largest clones (in terms of productive sequences), where n is:
#n = 100

# Processing all clones now...

# As long as they have at least a number min_clone_size of productive sequences, where:
min_clone_size = 20

dataset_ids = ['01107PB','01207PB','SFAPB']

for dataset in dataset_ids:
    partis_file_name = '../results/' + dataset + '_partition-cluster-annotations.csv'

    # ----- This section constrains the analysis to the largest n clones
    # # List of clone sizes
    # clone_size = []
    #
    # # Compute clone sizes
    # with open(partis_file_name, 'r') as partis_file:
    #     partis_csv = csv.DictReader(partis_file)
    #
    #     for row in partis_csv:
    #
    #         in_frame = row['in_frames'].split(':')
    #         has_stop_codons = row['stops'].split(':')
    #
    #         n_productive_seqs = [i for i in range(len(in_frame)) if in_frame[i] == 'True' and has_stop_codons[i] == 'False']
    #         n_productive_seqs = len(n_productive_seqs)
    #
    #         clone_size.append(n_productive_seqs)
    #
    # # Find the n largest clone sizes (multiple occurrences of a given size are counted multiple times):
    # # E.g. 5,5 and 5 are the 3 largest elements in [5,5,5,2,3]
    # n_largest_clone_sizes = nlargest(n, clone_size)
    #
    # largest_clone_indices = []
    # # Find indices (rows after header, starting at 0) of n largest clones with at least min_clone_size productive sequences:
    # # For each size value in n_largest_clone_sizes
    # for size in list(set(n_largest_clone_sizes)):
    #     if size >= min_clone_size:
    #         # Get the index  all clones of that size:
    #         clones_indices_with_size = [i for i in range(len(clone_size)) if clone_size[i] == size]
    #
    #         # Add clone indices to the list of the largest clone indices, to the limit of n:
    #         for index in clones_indices_with_size:
    #             if len(largest_clone_indices) < n:
    #                 largest_clone_indices.append(index)
    #
    # # Check that no duplicate indices were included:
    # assert len(largest_clone_indices) == len(list(set(largest_clone_indices)))

    with open(partis_file_name, 'r') as partis_file:
        partis_csv = csv.DictReader(partis_file)

        row_counter = 0

        for row in partis_csv:

            #if row_counter in largest_clone_indices:

            # If row is for a clone of successfully annotated sequences (e.g. has v gene info):
            if row['v_gene'] != '':


                # Find all clone sequences that were productive:
                in_frame = row['in_frames'].split(':')
                has_stop_codons = row['stops'].split(':')

                # Get indices of productive sequences
                p_indices = [i for i in range(len(in_frame)) if in_frame[i] == 'True' and has_stop_codons[i] == 'False']

                # Number of productive sequences:
                n_productive = len(p_indices)

                # If there is at least one productive sequence:
                if n_productive >= 1:

                    # Find sequence identifiers
                    seq_ids = row['unique_ids'].split(':')
                    seq_ids = [seq_ids[i] for i in range(len(seq_ids)) if i in p_indices]

                    assert len(seq_ids) == n_productive

                    # Find sequences:
                    clone_sequences = row['input_seqs'].split(':')
                    clone_sequences = [clone_sequences[i] for i in range(len(clone_sequences)) if i in p_indices]

                    # Remove 'N' nucleotides introduced by Partis at beginning of seq.
                    for i in range(len(clone_sequences)):
                        clone_sequences[i] = re.sub(r'N*[^ATCG]', '', clone_sequences[i], count=1)
                        clone_sequences[i] = re.sub(r'N*[^ATCG]\Z', '', clone_sequences[i], count=1)

                        #N_string = re.search(r'N*[^ATCG]', clone_sequences[i])
                        #if N_string != None:
                        #    clone_sequences[i] = clone_sequences[i].replace(N_string.group(), '')

                    # Find naive sequence
                    naive_seq = row['naive_seq']

                    # Remove 'N' nucleotides introduced by Partis at beginning or end of sequence:
                    naive_seq = re.sub(r'N*[^ATCG]', '',naive_seq, count = 1)
                    naive_seq = re.sub(r'N*[^ATCG]\Z', '', naive_seq, count=1)

                    # Write sequences to clone-specific fasta file:
                    clone_file_name = '../results/clones/' + dataset + '_clone_' + str(row_counter) + '.fasta'
                    print 'clone_' + str(row_counter)

                    with open(clone_file_name, 'w') as clone_fasta_file:
                        # Write naive sequence
                        clone_fasta_file.write('>NAIVE\n')
                        clone_fasta_file.write(naive_seq + '\n\n')

                        # Write observed sequences:
                        for i in range(len(seq_ids)):
                            clone_fasta_file.write('>' + seq_ids[i] + '\n')
                            clone_fasta_file.write(clone_sequences[i] + '\n\n')

            row_counter += 1