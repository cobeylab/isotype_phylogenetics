import sys
from dendropy import Tree
import re
import csv

def main(argv):

    tree_file_path = str(argv[1])


    # Get clone id from tree file path:
    clone_id = re.search(r'[0-9 | A-Z]*_clone_[0-9]*', tree_file_path).group()

    # Get file path for isotype annotation:
    isotype_file_path = '../results/clones/' + clone_id + '_isotypes.csv'

    # Generate annotated tree file path
    annotated_tree_path = '../results/trees/'  + clone_id + '_annotated_tree.nex'


    # Read tree newick string as a tree object using dendropy:
    with open(tree_file_path, 'rU') as tree_file:
        tree_string = tree_file.next()
        #print tree_string

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

    with open(annotated_tree_path, 'w') as output_file:
        tree.write(output_file, schema='nexus')

    #tree.write(path = annotated_tree_path, schema="nexus")
    #Tree.write(tree, path = annotated_tree_path, schema = 'nexus')

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)