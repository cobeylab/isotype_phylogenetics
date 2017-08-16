#Read .nex files, get a list of isotypes, shuffle the isotypes list and reassign isotypes to nodes for "numPermut" times.
import random
from dendropy import Tree

#dataDir = '../results/clones/'
#treeDir = '../results/trees/'
dataDir = './'
treeDir = './'

clone_id = '01207PB_clone_2_annotated_tree.nex'

numPermut = 10

def get_tree(treefName, dataDir):
	tree = Tree.get(path = dataDir+treefName, schema='nexus')
	assert isinstance(tree, Tree)
	return tree
	
def get_isotypes(tree):
	isotypes = []
	for node in tree.leaf_nodes():
		isotype = node.annotations.get_value('isotype')
		if isotype != None:
			isotypes.append(isotype)
	return isotypes
	
def shuffle_and_reassign(tree, isotypes):
	random.shuffle(isotypes)
	idx = 0
	for node in tree.leaf_nodes():
		if node.taxon.label != 'NAIVE':
			node.annotations['isotype'] = isotypes[idx]
			idx += 1
	return tree
	

tree = get_tree(clone_id, treeDir)
#mean_subTree_size = calc_mean_subTree_size(tree)

isotypes = get_isotypes(tree)

for n in range(numPermut):
	tree = shuffle_and_reassign(tree, isotypes)
	print (tree.leaf_nodes()[10].annotations['isotype'])
	#mean_subTree_size = calc_mean_subTree_size(tree)


