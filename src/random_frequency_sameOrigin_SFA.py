#generate shuffled datasets to test frequency of clones having both IgG and IgA is as expeced
import csv
import re
import glob
import random

numPermut = 1000

dataDir = '../results/clones/'
outf = open(dataDir+"frequency_sameOrigin_random_SFA.csv", "w")
outf.write("dataset,IgGIgA_together,IgG_noIgA,IgA_noIgG\n")

def read_isotypes(isotype_csv):
	isotypes_clone = []	
	for sequence in isotype_csv:
		isotype = sequence['isotype'][0:3]
		isotypes_clone.append(isotype)
	
	return isotypes_clone

def shuffle_and_reassign(isoList, clone_lens):

	random.shuffle(isoList)

	isotypes_clones = []
	idx = 0
	for cdx in range(len(clone_lens)):
		clone_len = clone_lens[cdx]
		isotypes = []
		for i in range(clone_len):
			isotypes.append(isoList[idx])
			idx += 1
		isotypes_clones.append(isotypes)
		
	return isotypes_clones

	
datasets = ["01107"]

for dataset in datasets:

	clone_lens = []
	isotypes_ds = []
	
	for clone_isotypes in glob.glob(dataDir+dataset+"*isotypes.csv"):
		cloneID = re.search(r'[0-9 | A-Z]*_clone_[0-9]*', clone_isotypes).group()
		cloneID = cloneID.split("_clone_")[1]
		
		isotype_file = open(clone_isotypes, "r")
		isotype_csv = csv.DictReader(isotype_file)

		isotypes = read_isotypes(isotype_csv)
		clone_lens.append(len(isotypes))
		isotypes_ds += isotypes

	for n in range(numPermut):
		isotypes_clones = shuffle_and_reassign(isotypes_ds, clone_lens)
		total_clones = 0
		sameOrigin_clones = 0
		IgG_noIgA_clones = 0
		IgA_noIgG_clones = 0
		
		for isotypes in isotypes_clones:
		
			total_clones += 1
				
			if 'IgG' in isotypes and 'IgA' in isotypes:
				sameOrigin_clones += 1
			elif 'IgG' in isotypes and 'IgA' not in isotypes:
				IgG_noIgA_clones += 1
			elif 'IgA' in isotypes and 'IgG' not in isotypes:
				IgA_noIgG_clones += 1
		
		sameOrigin_freq = 1.0*sameOrigin_clones/total_clones
		IgG_noIgA_freq = 1.0*IgG_noIgA_clones/total_clones
		IgA_noIgG_freq = 1.0*IgA_noIgG_clones/total_clones
		
		printout = dataset + ","
		printout += str(sameOrigin_freq) +","
		printout += str(IgG_noIgA_freq) + ","
		printout += str(IgA_noIgG_freq) + "\n"
		outf.write(printout)

outf.close()
