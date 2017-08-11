#!/usr/bin/python
import sys
from dendropy import Tree
import re
import csv
from copy import deepcopy


# LOADING TEST TREE
with open('test_tree.nex', 'rU') as tree_file:
    tree = Tree.get_from_stream(tree_file, schema='nexus')


def max_isotype_subtree(tree, node):
    '''Returns the maximum subtree containing the query node (dendropy object) and only terminal nodes with same isotype
    '''

    node_isotype = node.annotations.get_value('isotype')

    # Initialize focal ancestor as query node
    focal_ancestor = node

    # True if all descendants of focal ancestor are the same isotype as the query node.
    all_same_isotype = True

    # Go up the chain of ancestors, stop as soon an ancestor with at least 1 descendant of a different isotype is found
    while all_same_isotype == True:

        # Find all descendants of focal ancestor (excluding internal nodes with filter_fn)
        descendants = focal_ancestor.ageorder_iter(filter_fn = lambda nd: nd.is_leaf())

        descendant_isotypes = [d.annotations.get_value('isotype') for d in descendants]

        discordant_isotypes = [isot for isot in descendant_isotypes if isot != node_isotype]

        # If all descendants of focal ancestor have the same isotype as the query node,
        if len(discordant_isotypes) == 0:

            # Make focal ancestor the maximum subtree ancestor so far:
            max_subtree_ancestor = focal_ancestor
            #print max_subtree_ancestor

            # Move up the tree, make focal ancestor's parent node the focal ancestor
            focal_ancestor = focal_ancestor.parent_node

        # Else, if focal ancestor's descendants include a discordant isotype:
        else:
            # Stop searching
            all_same_isotype = False

    # Extract subtree rooted at maximum subtree ancestor:
    max_node_descendants = max_subtree_ancestor.ageorder_iter()
    max_node_descendants = [nd for nd in max_node_descendants]

    node_filter_fn = lambda nd: nd in max_node_descendants
    output_tree = tree.extract_tree(node_filter_fn=node_filter_fn)

    # Copy annotations over to output tree:
    for nd in output_tree.leaf_nodes():
        original_node = tree.find_node_with_taxon_label(nd.taxon.label)

        nd.annotations.add_new('isotype', original_node.annotations.get_value('isotype'))

    # Assert output tree only has sequences with the same isotype as the query node
    output_tree_isotypes = [nd.annotations.get_value('isotype') for nd in output_tree.leaf_nodes()]

    assert output_tree_isotypes.count(node_isotype) == len(output_tree_isotypes)

    return output_tree

