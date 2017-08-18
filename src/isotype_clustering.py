#!/usr/bin/python
import sys
import os
from dendropy import Tree
import re
import csv
from copy import deepcopy
from numpy import mean


# Reading test tree:
with open(os.path.dirname(os.path.realpath(sys.argv[0])) + '/test_tree.nex', 'rU') as tree_file:
    test_tree = Tree.get_from_stream(tree_file, schema='nexus')

def max_isotype_subtree(tree, node):
    """Returns the maximum subtree containing the query node (dendropy object) and only terminal nodes with same isotype.
    The maximum isotype subtree is the largest tree that contains the node and only includes nodes with the same isotype as the focal node.

    The maximum isotype subtree for sequence H5ULPI203FRBWQ in the test tree contains only the sequence itself:

    >>> len(max_isotype_subtree(test_tree, test_tree.find_node_with_taxon_label('H5ULPI203FRBWQ')).nodes())
    1

    The same for sequence H5ULPI203FRBWQ:

    >>> len(max_isotype_subtree(test_tree, test_tree.find_node_with_taxon_label('H5ULPI203FMLK4')).nodes()) == 1
    1

    The maximum isotype subtree for sequence H5ULPI203GN1VH contains 269 sequences:

    >>> len(max_isotype_subtree(test_tree, test_tree.find_node_with_taxon_label('H5ULPI203GN1VH')).leaf_nodes())
    269
    """

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


def mean_isotype_subtree_size(tree, isotype_list = ['IgG','IgA']):
    """
    269 of the IgA sequences are in a single max. isotype subtree of size 269, 2 others are in subtrees of size 1
    This yields a mean subtree size of 267.0221402214022:

    >>> mean_isotype_subtree_size(test_tree, isotype_list = ['IgG','IgA'])['IgA'] == mean([269] * 269 + [1,1])
    True

    The IgG sequences of the test tree are distributed this way (by visual inspection):
    [16] * 16 + [4] * 4 + [2] * 2 + [4] * 4 + [2] * 2 + [1] * 2 + [1] + [2] * 2 + [2] * 2

    I.e. 16 sequences in a 16-sequence maximum subtree, 4 in a 4-sequence subtree, 2 in a 2-seq subtree, 4 in another 4-sequence subtree, etc.

    Yielding a mean subtree size of 8.7714285714285722:
    >>> mean_isotype_subtree_size(test_tree, isotype_list = ['IgG','IgA'])['IgG']
    8.7714285714285722
    """
    # Dictionary storing (terminal) node isotypes:
    isotype_dic = {}

    # Dictionary storing size of maximum subtree containing each node and only nodes of identical subtype
    isotype_subtree_size = {}

    # Retrieve isotype of all terminal nodes:
    for node in tree.leaf_nodes():
        isotype_dic[node.taxon.label] = node.annotations.get_value('isotype')

    # Analyze all nodes except NAIVE to find their maximum isotype subtree:
    nodes_to_analyze = [tree.find_node_with_taxon_label(seq_id) for seq_id in isotype_dic.keys() if seq_id != 'NAIVE']

    # Analyze only nodes with isotypes in isotype_list (ignores nodes with unknown isotype)
    nodes_to_analyze = [node for node in nodes_to_analyze if node.annotations.get_value('isotype') in isotype_list]

    while len(nodes_to_analyze) > 0:

        # Take first node in list of nodes to analyze
        focal_node = nodes_to_analyze[0]

        # Find its maximum isotype subtree with max_isotype_subtree function:
        max_subtree = max_isotype_subtree(tree, focal_node)

        subtree_size = len(max_subtree.leaf_nodes())

        isotype_subtree_size[focal_node.taxon.label] = subtree_size

        nodes_to_analyze.remove(focal_node)

        # List of seq ids contained in max subtree
        ids_in_subtree = [nd.taxon.label for nd in max_subtree.leaf_nodes()]

        # For any node contained in the focal node subtree, the subtree is also the max. subtree for that node
        nodes_in_focal_node_subtree = []

        for node in nodes_to_analyze:
            if node.taxon.label in ids_in_subtree:
                isotype_subtree_size[node.taxon.label] = subtree_size

                # Keep track of which nodes were incidentally analyzed due to being in the focal node subtree
                nodes_in_focal_node_subtree.append(node)

        # Remove nodes in the focal node subtree from list of nodes to analyze
        nodes_to_analyze = [node for node in nodes_to_analyze if node not in nodes_in_focal_node_subtree]

    # Compute average isotype subtree size; total and by isotype
    mean_subtree_size = {}
    mean_subtree_size['all'] = mean(isotype_subtree_size.values())

    for itype in isotype_list:
        values = [isotype_subtree_size[key] for key in isotype_subtree_size.keys() if isotype_dic[key] == itype]
        mean_subtree_size[itype] = mean(values)

    return mean_subtree_size

#if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

