#Read .nex files, get a list of isotypes, shuffle the isotypes list and reassign isotypes to nodes for "numPermut" times.
import random
import re
import sys
import numpy
import math
from dendropy import Tree
from isotype_clustering import mean_isotype_subtree_size
from isotype_clustering import test_tree

numPermut = 1000

def get_tree(treefName):
	tree = Tree.get_from_path(src=treefName, schema='nexus')
	return tree

def get_isotypes(tree):
	#get isotypes of all tips except for NAIVE sequence and return a list of isotypes
	isotypes = []
	for node in tree.leaf_nodes():
		isotype = node.annotations.get_value('isotype')
		if isotype != None:
			isotypes.append(isotype)
	return isotypes

def shuffle_and_reassign(tree, isotypes):
	#shuffle a list "isotypes" and reassign shuffled isotype annotations to tips except for NAIVE sequence
	random.shuffle(isotypes)
	idx = 0
	for node in tree.leaf_nodes():
		if node.taxon.label != 'NAIVE':
			node.annotations['isotype'] = isotypes[idx]
			idx += 1
	return tree

def get_pValue(mean_subtree_size_obs_isotype, isotype_subtree_rdm):
	#isotype_subtree_rdm is a list of means of isotype subtree sizes from randomizations
	p = 0
	for i in range(len(isotype_subtree_rdm)):
		if isotype_subtree_rdm[i] >= mean_subtree_size_obs_isotype:
			p += 1
	p = 1.0*p/1000
	return p

def main(treefName, outfName):

	tree = get_tree(treefName)

	isotypes = get_isotypes(tree)

	n_IgG = isotypes.count('IgG')
	n_IgA = isotypes.count('IgA')

	#test for clustering only if the clone has at least 2 of IgG and at least 2 IgA
	if n_IgG < 2 or n_IgA < 2:
		return

	outf = open(outfName, "w")

	#get mean isotype subtree size of the observed clone
	mean_subtree_size_obs = mean_isotype_subtree_size(tree)

	#get mean isotype subtree sizes of permutaitons
	IgG_subtree_rdm = []
	IgA_subtree_rdm = []
	all_subtree_rdm = []

	for n in range(numPermut):
		tree = shuffle_and_reassign(tree, isotypes)
		mean_subtree_size = mean_isotype_subtree_size(tree)
		IgG_subtree_rdm.append(mean_subtree_size['IgG'])
		IgA_subtree_rdm.append(mean_subtree_size['IgA'])
		all_subtree_rdm.append(mean_subtree_size['all'])

	#write an output file
	#for now, this contains both 1000 results of mean subtree sizes from 1000 randomizations and a summary line

	#write 1000 results from 1000 randomiazations
	outf.write('IgG_rdm,IgA_rdm,all_rdm\n')
	for n in range(numPermut):
		outf.write(str(IgG_subtree_rdm[n])+","+str(IgA_subtree_rdm[n])+","+str(all_subtree_rdm[n])+"\n")

	#write a summary line
	mean_IgG_subtree_rdm = numpy.mean(IgG_subtree_rdm)
	mean_IgA_subtree_rdm = numpy.mean(IgA_subtree_rdm)
	SD_IgG_subtree_rdm = numpy.std(IgG_subtree_rdm, ddof=1)
	SD_IgA_subtree_rdm = numpy.std(IgA_subtree_rdm, ddof=1)
	p_IgG_subtree_rdm = get_pValue(mean_subtree_size_obs['IgG'], IgG_subtree_rdm)
	p_IgA_subtree_rdm = get_pValue(mean_subtree_size_obs['IgA'], IgA_subtree_rdm)

	outf.write("clone_id,n_IgG,n_IgA,mean_IgG_subtree_size,mean_IgA_subtree_size,mean_IgG_subtree_rdm,mean_IgA_subtree_rdm,SD_IgG_subtree_rdm,SD_IgA_subtree_rdm,p_IgG_subtree_rdm,p_IgA_subtree_rdm\n")
	clone_id = re.search(r'[0-9 | A-Z]*_clone_[0-9]*', treefName).group()
	oneline = clone_id + ","
	oneline += str(n_IgG) + "," + str(n_IgA) + ","
	oneline += str(mean_subtree_size_obs['IgG']) + "," + str(mean_subtree_size_obs['IgA']) + ","
	oneline += str(mean_IgG_subtree_rdm) + "," + str(mean_IgA_subtree_rdm) + ","
	oneline += str(SD_IgG_subtree_rdm) + "," + str(SD_IgA_subtree_rdm) + ","
	oneline += str(p_IgG_subtree_rdm) + "," + str(p_IgA_subtree_rdm) + "\n"
	outf.write(oneline)
	outf.close()

if __name__ == '__main__':
	treefName = sys.argv[1]
	outfName = sys.argv[2]
	main(treefName, outfName)
