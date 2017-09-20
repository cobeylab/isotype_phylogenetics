#Sizes of clones having only IgG or having only IgA 
#Is there difference in clone-size distributions between IgG only clones and IgA only clones?

import csv
import re
import glob

dataDir = '../results/clones/'

def count_isotypes(isotype_csv):
	
	numIgG = 0
	numIgA = 0
	cloneSize = 0	
	for sequence in isotype_csv:
		cloneSize += 1
		isotype = sequence['isotype'][0:3]
		if isotype == 'IgG':
			numIgG += 1
		elif isotype == 'IgA':
			numIgA += 1
		
	effSize = numIgG + numIgA
	
	return numIgG, numIgA, effSize, cloneSize

outf = open(dataDir+"frac_clone_with_effSize.txt", "w")
outf.write("dataset,effSize,onlyIgA,onlyIgG,onlyIgA_or_onlyIgG\n")
datasets = ["01107", "01207", "SFA"]

for dataset in datasets:
	IgA_clone_sizes = []
	IgA_num_clones = []
	IgG_clone_sizes = []
	IgG_num_clones = []
	IgGIgA_num_clones = []

	for clone_isotypes in glob.glob(dataDir+dataset+"*isotypes.csv"):
		cloneID = re.search(r'[0-9 | A-Z]*_clone_[0-9]*', clone_isotypes).group()
		cloneID = cloneID.split("_clone_")[1]
		
		isotype_file = open(clone_isotypes, "r")
		isotype_csv = csv.DictReader(isotype_file)

		numIgG, numIgA, effSize, cloneSize = count_isotypes(isotype_csv)	
		
		if effSize == 0:
			continue

		if numIgG == 0:
			IgA_clone_sizes.append(effSize)
		elif numIgA == 0:
			IgG_clone_sizes.append(effSize)

	temp_IgA = IgA_clone_sizes[:]
	temp_IgG = IgG_clone_sizes[:]
	temp_IgGIgA = IgG_clone_sizes + IgA_clone_sizes			

	clone_sizes = sorted(list(set(IgG_clone_sizes + IgA_clone_sizes)))

	for size in clone_sizes:
		num = temp_IgA.count(size)
		IgA_num_clones.append(1.0*num/len(temp_IgA))
	
	for size in clone_sizes:
		num = temp_IgG.count(size)
		IgG_num_clones.append(1.0*num/len(temp_IgG))

	for size in clone_sizes:
		num = temp_IgGIgA.count(size)
		IgGIgA_num_clones.append(1.0*num/len(temp_IgGIgA))
	
	for i in range(len(clone_sizes)):
		oneline = dataset + ","
		oneline += str(clone_sizes[i])+","+str(IgG_num_clones[i])+","+str(IgA_num_clones[i])+","
		oneline += str(IgGIgA_num_clones[i])
		outf.write(oneline+"\n")	

outf.close()	
	
		
