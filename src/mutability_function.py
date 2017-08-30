#!/usr/bin/python

"""Implements a function for computing several mutability metrics for a DNA sequence based on a set of specified sequence partitions
"""

__author__ = "Marcos Vieira (mvieira@uchicago.edu)"

# imports

import csv
import os
import itertools
import numpy as np

# Get location of this script (allows files to be read by this script when calling this script from other scripts)
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# constants
# Test sequence:

test_sequence = "caggttcagctggtgcagtctggagct---gaggtgaagaagcctggggcctcagtgaaggtctcctgcaaggcttctggttacaccttt------------accagctatggtatcagctgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcagcgcttac------aatggtaacacaaactatgcacagaagctccag---ggcagagtcaccatgaccacagacacatccacgagcacagcctacatggagctgaggagcctgagatctgacgacacggccgtgtattactgtgcgagaga"
test_partition = [3, 29, 95]

test_WRCH = 'TGCTCCCCTGCT'
test_DGYW = 'GGTTGGGGAGCT'
test_WA = 'TAACCCTA'

test_sequence = test_sequence.upper()

# mutability scores from the S5F model by Yaari et al:
csvfile = open(__location__ + '/Yaari_fivemer_scores.csv', 'r')
csvread = csv.reader(csvfile, delimiter=' ', quotechar='|')
S5F = {}
firstline = True
for row in csvread:
    if firstline:
        firstline = False
    else:
        csvrow = row
        motif = csvrow[0].replace('"', '')
        mutability = csvrow[1].replace('"', '')
        S5F[motif] = float(mutability)

# constants:
W = {'A', 'T'}
Y = {'C', 'T'}
R = {'G', 'A'}
S = {'G', 'C'}
D = {'A', 'G', 'T'}
H = {'T', 'A', 'C'}
K = {'T', 'G'}
M = {'C', 'A'}
V = {'A', 'C', 'G'}
B = {'C', 'G', 'T'}
N = {'A', 'T', 'C', 'G'}
X = {'A', 'T', 'C', 'G'}

ambiguity_code = {'W': W, 'Y': Y, 'R': R, 'S': S,
                  'D': D, 'H': H, 'K': K, 'M': M,
                  'V': V, 'B': B, 'N': N, 'X': X}

# mutability scores from the background indep. component of the  7mer model of Elhanati et al:
csvfile = open(__location__ + '/Elhanati_files/7mer_freq_indep.csv', 'r')
csvread = csv.reader(csvfile, delimiter=',')
SevenMer = {}
for row in csvread:
    csvrow = row
    motif = csvrow[0]
    mutability = csvrow[1]
    SevenMer[motif] = float(mutability)


# Function for computing mean S5F score:

def compute_mean_S5F(seq):
    mean_S5F = 0

    # Set of ambiguities:
    ambiguities = {'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N'}

    # For each nucleotide, except the first two and the last two:
    for i in range(2, (len(seq) - 2)):
        motif = seq[i - 2:i + 3]

        # Check if motif contains ambiguous nucleotides
        # If it doesn't, find S5F from dictionary:
        if len(ambiguities.intersection(set(motif))) == 0:
            mean_S5F = mean_S5F + S5F[motif]

        else:
            # If it does, find the ambiguous positions:
            positions = [p for p in range(len(motif)) if len(ambiguities.intersection(set(motif[p]))) > 0]

            # Find the exact ambiguity at each site:
            ambiguous_nts = [motif[p] for p in positions]

            # Find the sets of possible nucleotides for each site:
            alternative_nt_sets = [ambiguity_code[nt] for nt in ambiguous_nts]

            # Find all possible choices of alternative nucleotides for each position:
            replacements = list(itertools.product(*alternative_nt_sets))

            set_lengths = [len(nt_set) for nt_set in alternative_nt_sets]
            n_replacements = np.prod(np.array(set_lengths))

            assert len(replacements) == n_replacements

            S5F_values = []
            for replacement in replacements:
                possible_motif = list(motif)
                for i in range(len(replacement)):
                    possible_motif[positions[i]] = replacement[i]
                possible_motif = ''.join(possible_motif)
                S5F_values.append(S5F[possible_motif])

            mean_value = float(sum(S5F_values)) / len(S5F_values)
            mean_S5F = mean_S5F + mean_value

    mean_S5F = mean_S5F / float(len(seq) - 4)
    return mean_S5F


# Function for computing mean 7mer score:

def compute_mean_7M(seq):
    # Set of ambiguities:
    ambiguities = {'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N'}

    mean_7M = 0

    # For each nucleotide, except the first three and the last three:
    for i in range(3, (len(seq) - 3)):
        motif = seq[i - 3:i + 4]

        # Check if motif contains ambiguous nucleotides
        # If it doesn't, find 7-mer mutability from dictionary:
        if len(ambiguities.intersection(set(motif))) == 0:
            mean_7M = mean_7M + SevenMer[motif]
        else:
            # If it does, find the ambiguous positions:
            positions = [p for p in range(len(motif)) if len(ambiguities.intersection(set(motif[p]))) > 0]

            # Find the exact ambiguity at each site:
            ambiguous_nts = [motif[p] for p in positions]

            # Find the sets of possible nucleotides for each site:
            alternative_nt_sets = [ambiguity_code[nt] for nt in ambiguous_nts]

            # Find all possible choices of alternative nucleotides for each position:
            replacements = list(itertools.product(*alternative_nt_sets))

            set_lengths = [len(nt_set) for nt_set in alternative_nt_sets]
            n_replacements = np.prod(np.array(set_lengths))

            assert len(replacements) == n_replacements

            SevenMer_values = []
            for replacement in replacements:
                possible_motif = list(motif)
                for i in range(len(replacement)):
                    possible_motif[positions[i]] = replacement[i]
                possible_motif = ''.join(possible_motif)
                SevenMer_values.append(SevenMer[possible_motif])

            mean_value = float(sum(SevenMer_values)) / len(SevenMer_values)
            mean_7M = mean_7M + mean_value

    mean_7M = mean_7M / float(len(seq) - 6)
    return mean_7M


# Function for -counting- 7mers with non-zero mutability
def count_nonzero_7M(seq):
    # Set of ambiguities:
    ambiguities = {'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N'}

    nonzero_7M = 0
    # For each nucleotide, except the first three and the last three:
    for i in range(3, (len(seq) - 3)):
        # If motif is non-ambiguous
        if len(ambiguities.intersection(set(seq[i - 3:i + 4]))) == 0:
            if SevenMer[seq[i - 3:i + 4]] > 0:
                nonzero_7M = nonzero_7M + 1
    return nonzero_7M


# Functions for counting AID hotspots:


def count_WRCH(seq):
    n_WRCH = 0
    # For each nucleotide, except the first two and the last one...
    for i in range(2, (len(seq) - 1)):
        if (seq[i - 2] in W) and (seq[i - 1] in R) and (seq[i] == 'C') and (seq[i + 1] in H):
            n_WRCH += 1
    return n_WRCH


def count_DGYW(seq):
    n_DGYW = 0
    # For each nucleotide, except the first and the last two...
    for i in range(1, (len(seq) - 2)):
        if (seq[i - 1] in D) and (seq[i] == 'G') and (seq[i + 1] in Y) and (seq[i + 2] in W):
            n_DGYW += 1
    return n_DGYW


# Functions for counting POL hotspots:


def count_WA(seq):
    n_WA = 0
    # For each nucleotide, except the first one:
    for i in range(1, len(seq)):
        if (seq[i - 1] in W) and (seq[i] == 'A'):
            n_WA += 1
    return n_WA


def count_TW(seq):
    n_TA = 0
    # For each nucleotide, except the last one:
    for i in range(0, (len(seq) - 1)):
        if (seq[i] == 'T') and (seq[i + 1] in W):
            n_TA += 1
    return n_TA


# Functions for counting AID coldspots:
def count_SYC(seq):
    n_SYC = 0
    # For each nucleotide, except the first two
    for i in range(2, len(seq)):
        if (seq[i - 2] in S) and (seq[i - 1] in Y) and (seq[i] == "C"):
            n_SYC += 1
    return n_SYC


def count_GRS(seq):
    n_GRS = 0
    # For each nucleotide, except the last two:
    for i in range(0, (len(seq) - 2)):
        if (seq[i] == "G") and (seq[i + 1] in R) and (seq[i + 2] in S):
            n_GRS += 1
    return n_GRS


def count_trinucleotide_CS_coding(seq):
    n_trinucleotide_CS_coding = 0

    # Coding strand coldspot trinucleotides
    coding_CS = {"TTC", "CAC", "GGC", "GAC"}

    # For each nucleotide, except the first two:
    for i in range(2, len(seq)):
        tri = seq[i - 2:i + 1]
        if (tri in coding_CS):
            n_trinucleotide_CS_coding += 1
    return n_trinucleotide_CS_coding


def count_trinucleotide_CS_noncoding(seq):
    n_trinucleotide_CS_noncoding = 0
    # non-coding strand coldspot trinucleotides
    noncoding_CS = {"GAA", "GTG", "GCC", "GTC"}
    # For each nucleotide, except the last two:
    for i in range(0, (len(seq) - 2)):
        tri = seq[i:i + 3]
        if (tri in noncoding_CS):
            n_trinucleotide_CS_noncoding += 1
    return n_trinucleotide_CS_noncoding


# Function for computing overlapping hotspots (Wei et al. 2015 PNAS):
# The brute-force approach was just inertia, could use regex.
def count_OHS(seq):
    n_OHS = 0
    # For each nucleotide, except the first two and the last one:
    for i in range(2, (len(seq) - 1)):
        if (seq[i - 2] + seq[i - 1] + seq[i] + seq[i + 1]) == 'AGCT':
            n_OHS += 1
    return n_OHS


# Main function:

def seq_mutability(seq, partition_points=None):
    """Calculates several mutability metrics for a specified set of contiguous partitions in a DNA sequence (seq). The partitions are specified by argument partition_points. Each element of the list partition_points specifies the position in the sequence (NOT the Python index) where a partition begins. The last element indicates the last position of the last partition (in case the user does not want to cover the entire sequence). If partition_points is set to 'None', the whole sequence is treated as a single partition. Deals with gaps.
    """
    seq = seq.upper()
    if partition_points is None:
        partition_points = [1, len(seq)]

    # Convert partition_points into Python indices (i.e. subtract 1):
    partition_points = [(partition_points[i] - 1)
                        for i in range(len(partition_points))]

    # Output list of lists.First internal list gives partition lengths.
    # Each subsequent list gives mutability scores for a partition
    results_lists = [[]]

    # From partition points, store ungapped partitions in a list:
    ungapped_partitions = []

    for i in (range(len(partition_points) - 1)):
        if i < (len(partition_points) - 2):
            partition = seq[partition_points[i]: partition_points[i + 1]]
        else:
            # For the last partition, the correct parsing expression is
            partition = seq[partition_points[i]: partition_points[i + 1] + 1]
            # Since the last point is the last position of the las partition.
        ungapped_partitions.append(partition.replace("-", ""))
        results_lists[0].append(len(partition.replace("-", "")))

    # Find and ungap right and left "margins" - regions outside partitions
    # If 1st p. point not at the beggining of seq, then there's a left margin
    if partition_points[0] > 0:
        left_margin = seq[0:partition_points[0]]
        left_margin = left_margin.replace("-", "")
    else:
        left_margin = ""

    # If last p. point not the last position of seq, there's a right margin:
    if partition_points[len(partition_points) - 1] < (len(seq) - 1):
        right_margin = seq[(partition_points[len(partition_points) - 1]) + 1:(len(seq) - 1)]
        right_margin = right_margin.replace("-", "")
    else:
        right_margin = ""

    # For each (ungapped) partition:
    for i in range(len(ungapped_partitions)):
        partition = ungapped_partitions[i]

        # Find neighbor nucleotides
        # If this is the leftmost partition
        if i == 0:
            if len(left_margin) > 0:
                l_neighbors = ['', '', left_margin[-1]]
            if len(left_margin) > 1:
                l_neighbors = ['', left_margin[-2], left_margin[-1]]
            if len(left_margin) > 2:
                l_neighbors = [left_margin[-3], left_margin[-2], left_margin[-1]]
            if len(left_margin) == 0:
                l_neighbors = ['', '', '']
            # If there are other partitions:
            if len(ungapped_partitions) > 1:
                r_neighbors = [ungapped_partitions[i + 1][0], ungapped_partitions[i + 1][1],
                               ungapped_partitions[i + 1][2]]
            else:
                r_neighbors = ['', '', '']

        # If this is the rightmost partition:
        elif i == (len(ungapped_partitions) - 1):
            if len(right_margin) > 0:
                r_neighbors = [right_margin[0], '', '']
            if len(right_margin) > 1:
                r_neighbors = [right_margin[0], right_margin[1], '']
            if len(right_margin) > 2:
                r_neighbors = [right_margin[0], right_margin[1], right_margin[2]]
            if len(right_margin) == 0:
                r_neighbors = ['', '', '']
            # If there are other partitions:
            if len(ungapped_partitions) > 1:
                l_neighbors = [ungapped_partitions[i - 1][-3], ungapped_partitions[i - 1][-2],
                               ungapped_partitions[i - 1][-1]]
            else:
                l_neighbors = ['', '', '']
        # If this partition has neighbor partitions on both sides:
        else:
            l_neighbors = [ungapped_partitions[i - 1][-3], ungapped_partitions[i - 1][-2],
                           ungapped_partitions[i - 1][-1]]
            r_neighbors = [ungapped_partitions[i + 1][0], ungapped_partitions[i + 1][1],
                           ungapped_partitions[i + 1][2]]

        # Calc. mutability metrics, appending neighbors where necessary:
        mean_S5F = compute_mean_S5F(l_neighbors[1] + l_neighbors[2] +
                                    partition +
                                    r_neighbors[0] + r_neighbors[1])

        mean_7M = compute_mean_7M(l_neighbors[0] + l_neighbors[1] + l_neighbors[2] +
                                  partition +
                                  r_neighbors[0] + r_neighbors[1] + r_neighbors[2])

        n_nonzero_7M = count_nonzero_7M(l_neighbors[0] + l_neighbors[1] + l_neighbors[2] +
                                        partition +
                                        r_neighbors[0] + r_neighbors[1] + r_neighbors[2])

        n_WRCH = count_WRCH(l_neighbors[1] + l_neighbors[2] +
                            partition + r_neighbors[0])

        n_DGYW = count_DGYW(l_neighbors[2] + partition +
                            r_neighbors[0] + r_neighbors[1])

        n_AIDHS = n_WRCH + n_DGYW

        n_WA = count_WA(l_neighbors[2] + partition)

        n_TW = count_TW(partition + r_neighbors[0])

        n_POLHS = n_WA + n_TW

        C_fraction = partition.count('C')
        C_fraction = C_fraction / float(len(partition))

        n_SYC = count_SYC(l_neighbors[1] + l_neighbors[2] + partition)
        n_GRS = count_GRS(partition + r_neighbors[0] + r_neighbors[1])
        trinuc_CS_coding = count_trinucleotide_CS_coding(l_neighbors[1] +
                                                         l_neighbors[2] +
                                                         partition)

        trinuc_CS_noncoding = count_trinucleotide_CS_noncoding(partition +
                                                               r_neighbors[0] +
                                                               r_neighbors[1])

        n_AIDCS = n_SYC + n_GRS + trinuc_CS_coding + trinuc_CS_noncoding

        n_OHS = count_OHS(l_neighbors[1] + l_neighbors[2] +
                          partition + r_neighbors[0])

        partition_results = [mean_S5F, mean_7M, n_nonzero_7M, n_AIDHS, n_POLHS, C_fraction, n_AIDCS, n_OHS]
        partition_results = {'mean_S5F': mean_S5F, 'mean_7M': mean_7M, 'n_nonzero_7M': n_nonzero_7M,
                             'HS': n_AIDHS, 'POLHS': n_POLHS, 'C_fraction': C_fraction, 'CS': n_AIDCS,
                             'OHS': n_OHS}

        results_lists.append(partition_results)
    return results_lists


def aggregated_mutability(seq, partition_points):
    '''Calls seq_mutability and aggregates resultsinto FRs (partitions 1,3,5) vs CDRs (partitions 2,4,6)
       Hotspots, coldspots, polymerase hotspots, overlapping HS and non-zero 7Ms are summed, others are averaged weighted by the length of each partition"
    '''
    mutability = seq_mutability(seq, partition_points)

    length_FR1, length_CDR1, length_FR2, length_CDR2, length_FR3, length_CDR3 = mutability[0]

    weight_FR1 = length_FR1 / float(length_FR1 + length_FR2 + length_FR3)
    weight_FR2 = length_FR2 / float(length_FR1 + length_FR2 + length_FR3)
    weight_FR3 = length_FR3 / float(length_FR1 + length_FR2 + length_FR3)

    weight_CDR1 = length_CDR1 / float(length_CDR1 + length_CDR2 + length_CDR3)
    weight_CDR2 = length_CDR2 / float(length_CDR1 + length_CDR2 + length_CDR3)
    weight_CDR3 = length_CDR3 / float(length_CDR1 + length_CDR2 + length_CDR3)

    combined_FR_S5F = mutability[1]['mean_S5F'] * weight_FR1 + mutability[3]['mean_S5F'] * weight_FR2 + mutability[5][
                                                                                                            'mean_S5F'] * weight_FR3
    combined_FR_7M = mutability[1]['mean_7M'] * weight_FR1 + mutability[3]['mean_7M'] * weight_FR2 + mutability[5][
                                                                                                         'mean_7M'] * weight_FR3
    combined_FR_nonzero_7M = mutability[1]['n_nonzero_7M'] + mutability[3]['n_nonzero_7M'] + mutability[5][
        'n_nonzero_7M']
    combined_FR_AIDHS = mutability[1]['HS'] + mutability[3]['HS'] + mutability[5]['HS']
    combined_FR_POLHS = mutability[1]['POLHS'] + mutability[3]['POLHS'] + mutability[5]['POLHS']
    combined_FR_C_fraction = mutability[1]['C_fraction'] * weight_FR1 + mutability[3]['C_fraction'] * weight_FR2 + \
                             mutability[5]['C_fraction'] * weight_FR3
    combined_FR_AIDCS = mutability[1]['CS'] + mutability[3]['CS'] + mutability[5]['CS']
    combined_FR_OHS = mutability[1]['OHS'] + mutability[3]['OHS'] + mutability[5]['OHS']

    combined_CDR_S5F = mutability[2]['mean_S5F'] * weight_CDR1 + mutability[4]['mean_S5F'] * weight_CDR2 + \
                       mutability[6]['mean_S5F'] * weight_CDR3
    combined_CDR_7M = mutability[2]['mean_7M'] * weight_CDR1 + mutability[4]['mean_7M'] * weight_CDR2 + mutability[6][
                                                                                                            'mean_7M'] * weight_CDR3
    combined_CDR_nonzero_7M = mutability[2]['n_nonzero_7M'] + mutability[4]['n_nonzero_7M'] + mutability[6][
        'n_nonzero_7M']
    combined_CDR_AIDHS = mutability[2]['HS'] + mutability[4]['HS'] + mutability[6]['HS']
    combined_CDR_POLHS = mutability[2]['POLHS'] + mutability[4]['POLHS'] + mutability[6]['POLHS']
    combined_CDR_C_fraction = mutability[2]['C_fraction'] * weight_CDR1 + mutability[4]['C_fraction'] * weight_CDR2 + \
                              mutability[6]['C_fraction'] * weight_CDR3
    combined_CDR_AIDCS = mutability[2]['CS'] + mutability[4]['CS'] + mutability[6]['CS']
    combined_CDR_OHS = mutability[2]['OHS'] + mutability[4]['OHS'] + mutability[6]['OHS']

    return {'FR_mutability': {'mean_S5F': combined_FR_S5F, 'mean_7M': combined_FR_7M,
                              'n_nonzero_7M': combined_FR_nonzero_7M,
                              'HS': combined_FR_AIDHS, 'POLHS': combined_FR_POLHS, 'C_fraction': combined_FR_C_fraction,
                              'CS': combined_FR_AIDCS, 'OHS': combined_FR_OHS},
            'CDR_mutability': {'mean_S5F': combined_CDR_S5F, 'mean_7M': combined_CDR_7M,
                               'n_nonzero_7M': combined_CDR_nonzero_7M,
                               'HS': combined_CDR_AIDHS, 'POLHS': combined_CDR_POLHS,
                               'C_fraction': combined_CDR_C_fraction,
                               'CS': combined_CDR_AIDCS, 'OHS': combined_CDR_OHS}
            }
