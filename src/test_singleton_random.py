import csv
import re
import glob
#dataDir = '../../../mvieira/isotype_phylogenetics/results/clones/'
dataDir = '../results/clones/'

outf = open(dataDir+"isotypeFrac_cloneSize.csv", "w")
singlef = open(dataDir+"singleton_fractions.csv", "w")

outf.write("dataset,cloneID,numIgG,numIgA,effSize,cloneSize,effFracIgG,effFracIgA,fracIgG,fracIgA\n")
singlef.write("dataset, numIgG, numIgA, fracIgG, fracIgA, numSingleG, numSingleA, fracSingleG, fracSingleA\n")

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

datasets = ["01107", "01207", "SFA"]

numIgG_all = 0
numIgA_all = 0
numSingleG_all = 0
numSingleA_all = 0

for dataset in datasets:
	numIgG_ds = 0 
	numIgA_ds = 0
	numSingleG_ds = 0 #number of IgG singleton
	numSingleA_ds = 0 #number of IgA singleton
	for clone_isotypes in glob.glob(dataDir+dataset+"*isotypes.csv"):
		cloneID = re.search(r'[0-9 | A-Z]*_clone_[0-9]*', clone_isotypes).group()
		cloneID = cloneID.split("_clone_")[1]
		
		isotype_file = open(clone_isotypes, "r")
		isotype_csv = csv.DictReader(isotype_file)

		numIgG, numIgA, effSize, cloneSize = count_isotypes(isotype_csv)	

		numIgG_ds += numIgG
		numIgA_ds += numIgA		
		if cloneSize == 1 and numIgG == 1:
			numSingleG_ds += 1
		elif cloneSize == 1 and numIgA == 1:
			numSingleA_ds += 1

		try:
			effFracIgG = 1.0*numIgG/effSize
		except ZeroDivisionError:
			effFracIgG = 'NA'
		try:
			effFracIgA = 1.0*numIgA/effSize
		except ZeroDivisionError:
			effFracIgA = 'NA'
			
		fracIgG = 1.0*numIgG/cloneSize
		fracIgA = 1.0*numIgA/cloneSize
		
		#write file for correlation test
		oneline = dataset+","+cloneID+","
		oneline += str(numIgG)+","+str(numIgA)+","+str(effSize)+","+str(cloneSize)+","
		oneline += str(effFracIgG)+","+str(effFracIgA)+","+str(fracIgG)+","+str(fracIgA)
		outf.write(oneline+"\n")
		
		isotype_file.close()
	
	#write singleton fraction file
	fracIgG_ds = 1.0*numIgG_ds/(numIgG_ds+numIgA_ds)
	fracIgA_ds = 1.0*numIgA_ds/(numIgG_ds+numIgA_ds)
	fracSingleG_ds = 1.0*numSingleG_ds/(numSingleG_ds+numSingleA_ds)
	fracSingleA_ds = 1.0*numSingleA_ds/(numSingleG_ds+numSingleA_ds)
	oneline = dataset+","+str(numIgG_ds)+","+str(numIgA_ds)+","+str(fracIgG_ds)+","+str(fracIgA_ds)+","
	oneline += str(numSingleG_ds)+","+str(numSingleA_ds)+","+str(fracSingleG_ds)+","+str(fracSingleA_ds)
	singlef.write(oneline+'\n')
	
	numIgG_all += numIgG_ds
	numIgA_all += numIgA_ds
	numSingleG_all += numSingleG_ds
	numSingleA_all += numSingleA_ds
	
outf.close()		

fracIgG_all = 1.0*numIgG_all/(numIgG_all + numIgA_all)
fracIgA_all = 1.0*numIgA_all/(numIgG_all + numIgA_all)
fracSingleG_all = 1.0*numSingleG_all/(numSingleG_all + numSingleA_all)
fracSingleA_all = 1.0*numSingleA_all/(numSingleG_all + numSingleA_all)
oneline = "all,"+str(numIgG_all)+","+str(numIgA_all)+","+str(fracIgG_all)+","+str(fracIgA_all)+","
oneline += str(numSingleG_all)+","+str(numSingleA_all)+","+str(fracSingleG_all)+","+str(fracSingleA_all)
singlef.write(oneline+"\n")

singlef.close()
		
