#!/usr/bin/python
import sys
from dendropy import Tree
import re
import csv

tree_file_path = '../results/trees/RAxML_bestTree.01207PB_clone_2'


# Get clone id from tree file path:
clone_id = re.search(r'[0-9 | A-Z]*_clone_[0-9]*', tree_file_path).group()

# Get file path for isotype annotation:
isotype_file_path = '../results/clones/' + clone_id + '_isotypes.csv'

# Generate annotated tree file path
annotated_tree_path = '../results/trees/'  + clone_id + '_annotated_tree.nex'


# Read tree newick string as a tree object using dendropy:
with open(tree_file_path, 'rU') as tree_file:
    tree_string = tree_file.next()

tree_string = tree_string.replace('\n','')

tree = Tree.get_from_string(tree_string, schema='newick')

# Read isotype csv file as a dictionary:
isotype = {}
with open(isotype_file_path, 'rU') as isotype_file:
    isotype_csv = csv.DictReader(isotype_file)
    for row in isotype_csv:
        # Remove space from ' isotype' after Kangchon changes script
        isotype[row['id']] = row['isotype'].translate(None, '1234')

# Annotate tree with isotype info:
for node in tree.leaf_nodes():
    seq_id = node.taxon.label

    # If sequence id is in the isotype dictionary (i.e. is not the NAIVE sequence or a short sequence w/ no isotype):
    if seq_id in isotype.keys():
        node.annotations.add_new('isotype', isotype[seq_id])

# Find MRCA of all IgG isotypes and MRCA of all IgA isotypes


# Extract subtree with MRCA

MRCA_descendants = MRCA.ageorder_iter()
MRCA_descendants = [nd for nd in MRCA_descendants]

node_filter_fn = lambda nd: nd in MRCA_descendants
tree2 = tree.extract_tree(node_filter_fn=node_filter_fn)




tree.write(path = annotated_tree_path, schema="nexus")










