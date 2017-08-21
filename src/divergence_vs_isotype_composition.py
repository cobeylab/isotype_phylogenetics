#!/usr/bin/python
"""Takes an annotated clone tree file, finds the mean and the maximum divergence of tips from the ancestor and  the IgA/IgG content"""
import sys
import os
from dendropy import Tree
import re
import csv
from copy import deepcopy
from numpy import mean

def main(argv):
    tree_file_path = str(argv[1])

    # Get clone id from tree file path:
    clone_id = re.search(r'[0-9 | A-Z]*_clone_[0-9]*', tree_file_path).group()

    # Get dataset id:
    #dataset_id = re.search(r'[^_]*', clone_id).group()

    # Generate temporary file with results for this clone
    temp_output_file_path = '../results/'  + clone_id + '_divergence_vs_composition_temp.csv'


    # Read tree
    tree = Tree.get_from_path(tree_file_path, 'nexus')

    # Naive seq. dendropy taxon identifier
    naive_taxon = [node.taxon for node in tree.leaf_nodes() if node.taxon.label.find('NAIVE') > -1][0]

    # Compute divergence of all sequences from naive sequence
    pmatrix = tree.phylogenetic_distance_matrix()

    divergs = [pmatrix.distance(node.taxon, naive_taxon) for node in tree.leaf_nodes() if node.taxon.label != 'NAIVE']


    max_divergence = str(max(divergs))
    mean_divergence = str(mean(divergs))
    cumulative_divergence = str(sum(divergs))


    tip_isotypes = [node.annotations.get_value('isotype') for node in tree.leaf_nodes()]

    n_igg = tip_isotypes.count('IgG')
    n_iga = tip_isotypes.count('IgA')

    # Naive sequence should not count towards number of sequences:
    n_sequences = len(tree.leaf_nodes()) - 1

    fraction_igg = str(float(n_igg) / n_sequences)
    fraction_iga = str(float(n_iga) / n_sequences)

    n_igg = str(n_igg)
    n_iga = str(n_iga)
    n_sequences = str(n_sequences)

    results_list = [clone_id, n_sequences, n_igg, n_iga, fraction_igg, fraction_iga,
                    mean_divergence, max_divergence, cumulative_divergence]

    with open(temp_output_file_path, 'w') as temp_output_file:
        temp_output_file.write('clone_id,n_sequences,n_igg,n_iga,fraction_igg,fraction_iga,')
        temp_output_file.write('mean_divergence,max_divergence,cumulative_divergence\n')

        temp_output_file.write(','.join(results_list))
        temp_output_file.write('\n')

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)

