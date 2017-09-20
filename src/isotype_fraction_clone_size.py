#Get fractions from clones having at least one IgG or IgA

#Outputs:
#	isotypeFrac_cloneSize.txt
#		Get the fraction of IgG and IgG within a clone along with its clone size
#		Larger clones have more/less IgG than expected?
#	mean_IgG_IgA_proportions_within_IgGIgAclones.csv
#	singleton_fractions.csv

import csv
import re
import glob

dataDir = '../results/clones/'

#write size and fraction of IgG and IgA in each clone
outf = open(dataDir+"isotypeFrac_cloneSize.txt", "w")
outf.write("dataset,cloneID,numIgG,numIgA,effSize,cloneSize,effFracIgG,effFracIgA,fracIgG,fracIgA\n")

#write mean proportion of IgG and IgA within clones having both IgG and IgA
propf = open(dataDir+"mean_IgG_IgA_proportions_within_IgGIgAclones.csv", "w")
propf.write("dataset,eff_IgG_prop,eff_IgA_prop,IgG_prop,IgA_prop\n")

#write fraction of singleton IgG and IgA
singlef = open(dataDir+"singleton_fractions.csv", "w")
singlef.write("dataset,numIgG,numIgA,fracIgG,fracIgA,numSingleG,numSingleA,fracSingleG,fracSingleA\n")


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

def write_singleton_fractions(dataset,numIgGs,numIgAs,cloneSizes):	
#singlef.write("dataset,numIgG,numIgA,fracIgG,fracIgA,numSingleG,numSingleA,fracSingleG,fracSingleA\n")
	sumIgG = sum(numIgGs)
	sumIgA = sum(numIgAs)
	fracIgG = 1.0*sumIgG/(sumIgG+sumIgA)
	fracIgA = 1.0*sumIgA/(sumIgG+sumIgA)
	numSingleG = 0
	numSingleA = 0
	for idx in range(len(numIgGs)):
		if (numIgGs[idx] == 1 and cloneSizes[idx] == 1) or (numIgAs[idx] == 1 and cloneSizes[idx] == 1): 
			numSingleG += numIgGs[idx]	
			numSingleA += numIgAs[idx]
	fracSingleG = 1.0*numSingleG/(numSingleG+numSingleA)
	fracSingleA = 1.0*numSingleA/(numSingleG+numSingleA)
	
	oneline = dataset + ","
	oneline += str(sumIgG) + "," + str(sumIgA) + "," 
	oneline += str(fracIgG) + "," + str(fracIgA) + ","
	oneline += str(numSingleG) + "," + str(numSingleA) + ","
	oneline += str(fracSingleG) + "," + str(fracSingleA)
	
	singlef.write(oneline+"\n")
	
	return oneline.split(",")
		

def write_prop_IgG_IgA_within_IgGIgAClone(dataset, effProps_IgG, effProps_IgA, props_IgG, props_IgA):

	effProps_IgG_in_GAclone = []
	effProps_IgA_in_GAclone = []
	props_IgG_in_GAclone = []
	props_IgA_in_GAclone = []
	
	for idx in range(len(effProps_IgG)):
		if effProps_IgG[idx] > 0 and effProps_IgA[idx] > 0:
			effProps_IgG_in_GAclone.append(effProps_IgG[idx])
			effProps_IgA_in_GAclone.append(effProps_IgA[idx])
			props_IgG_in_GAclone.append(props_IgG[idx])
			props_IgA_in_GAclone.append(props_IgA[idx])
	
	oneline = dataset+","
	oneline += str(sum(effProps_IgG_in_GAclone)/len(effProps_IgG_in_GAclone))+","
	oneline += str(sum(effProps_IgA_in_GAclone)/len(effProps_IgA_in_GAclone))+","
	oneline += str(sum(props_IgG_in_GAclone)/len(props_IgG_in_GAclone))+","
	oneline += str(sum(props_IgA_in_GAclone)/len(props_IgA_in_GAclone))
	propf.write(oneline+"\n")

	return oneline.split(",")


datasets = ["01107", "01207", "SFA"]
avg_prop = []
avg_single = []

for dataset in datasets:
	effProps_IgG = []
	effProps_IgA = []
	props_IgG = []
	props_IgA = []
	numIgGs = []
	numIgAs = []
	cloneSizes = []
	
	for clone_isotypes in glob.glob(dataDir+dataset+"*isotypes.csv"):
		cloneID = re.search(r'[0-9 | A-Z]*_clone_[0-9]*', clone_isotypes).group()
		cloneID = cloneID.split("_clone_")[1]
		
		isotype_file = open(clone_isotypes, "r")
		isotype_csv = csv.DictReader(isotype_file)

		numIgG, numIgA, effSize, cloneSize = count_isotypes(isotype_csv)	
		
		if effSize == 0:
			continue

		effFracIgG = 1.0*numIgG/effSize
		effFracIgA = 1.0*numIgA/effSize
		fracIgG = 1.0*numIgG/cloneSize
		fracIgA = 1.0*numIgA/cloneSize
		
		#write in outf file
		oneline = dataset+","+cloneID+","
		oneline += str(numIgG)+","+str(numIgA)+","+str(effSize)+","+str(cloneSize)+","
		oneline += str(effFracIgG)+","+str(effFracIgA)+","+str(fracIgG)+","+str(fracIgA)
		outf.write(oneline+"\n")
		
		#save proportions and sizes
		effProps_IgG.append(effFracIgG)
		effProps_IgA.append(effFracIgA)
		props_IgG.append(fracIgG)
		props_IgA.append(fracIgA)	
		numIgGs.append(numIgG)
		numIgAs.append(numIgA)
		cloneSizes.append(cloneSize)

		isotype_file.close()
		
	#write in propf file
	avg_prop.append(write_prop_IgG_IgA_within_IgGIgAClone(dataset, effProps_IgG, effProps_IgA, props_IgG, props_IgA))
	
	#write in singlef file
	avg_single.append(write_singleton_fractions(dataset,numIgGs,numIgAs, cloneSizes))

	
oneline = "average," 
for i in range(2, len(avg_prop[0])):
	avg = 0
	for j in range(len(avg_prop)):
		avg += float(avg_prop[j][i])
	avg = 1.0*avg/len(avg_prop)
	oneline += str(avg)
	if i<len(avg_prop[0])-1:
		oneline += ","
		
propf.write(oneline+"\n")

oneline = "average," 
for i in range(1, len(avg_single[0])):
	avg = 0
	for j in range(len(avg_single)):
		avg += float(avg_single[j][i])
	avg = 1.0*avg/len(avg_single)
	oneline += str(avg)
	if i<len(avg_single[0])-1:
		oneline += ","
		
singlef.write(oneline+"\n")

	
		
singlef.close()	
propf.close()
outf.close()		
