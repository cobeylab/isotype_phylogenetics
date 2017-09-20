import csv
import glob
import re
dataDir = '../../../mvieira/isotype_phylogenetics/results/clones/'

outf = open(dataDir+"clones_with_IgE.csv", "w")
outf.write("dataset,cloneID,numIgE,numIgG,numIgA,cloneSize\n")

meanf = open(dataDir+"IgE_summary.csv", "w")
meanf.write("dataset,frac_IgE_clones,frac_IgEIgG_clones,frac_IgEIgA_clones,frac_IgEIgGIgA_clones\n")

datasets = ["SFA"]

mean_eff_IgGs = []
mean_eff_IgAs = []
mean_IgGs = []
mean_IgAs = []

for dataset in datasets:
	print dataset
	#eff_IgG_props = []
	#eff_IgA_props = []
	#IgG_props = []
	#IgA_props = []
	numTotClones = 0
	numIgEClones = 0
	numIgEIgGClones = 0
	numIgEIgAClones = 0
	numIgEIgGIgAClones =0

	for clone_isotypes in glob.glob(dataDir+dataset+"*isotypes.csv"):
		
		isotype_file = open(clone_isotypes, "r")
		isotype_csv = csv.DictReader(isotype_file)
		cloneID = re.search(r'clone_[0-9]*', clone_isotypes).group()
		
		numTotClones += 1
		
		isotypes = []
		for sequence in isotype_csv:
			isotypes.append(sequence['isotype'][0:3])
						
		if 'IgE' in isotypes:
		
			numIgEClones += 1
		
			if 'IgG' in isotypes:
				numIgEIgGClones += 1
			if 'IgA' in isotypes:
				numIgEIgAClones += 1
			if 'IgG' in isotypes and 'IgA' in isotypes:
				numIgEIgGIgAClones += 1
	 
			numIgE = isotypes.count('IgE')
			numIgG = isotypes.count('IgG')
			numIgA = isotypes.count('IgA')
			#eff_IgG_prop = 1.0*numIgG/(numIgG+numIgA)
			#eff_IgA_prop = 1.0*numIgA/(numIgG+numIgA)
			#IgG_prop = 1.0*numIgG/len(isotypes)
			#IgA_prop = 1.0*numIgA/len(isotypes)

			#eff_IgG_props.append(eff_IgG_prop)
			#eff_IgA_props.append(eff_IgA_prop)
			#IgG_props.append(IgG_prop)
			#IgA_props.append(IgA_prop)
				
			printout = dataset + "," + cloneID + ","
			printout += str(numIgE) + "," + str(numIgG) + "," + str(numIgA) + "," + str(len(isotypes))
			outf.write(printout+"\n")
		
		isotype_file.close()
		
	#calculate mean for dataset
	#mean_eff_IgG = 1.0*sum(eff_IgG_props)/len(eff_IgG_props)
	#mean_eff_IgA = 1.0*sum(eff_IgA_props)/len(eff_IgA_props)
	#mean_IgG = 1.0*sum(IgG_props)/len(IgG_props)
	#mean_IgA = 1.0*sum(IgA_props)/len(IgA_props)
	
	#mean_eff_IgGs.append(mean_eff_IgG)
	#mean_eff_IgAs.append(mean_eff_IgA)
	#mean_IgGs.append(mean_IgG)
	#mean_IgAs.append(mean_IgA)	

	frac_IgE = 1.0*numIgEClones/numTotClones
	frac_IgEIgG = 1.0*numIgEIgGClones/numIgEClones
	frac_IgEIgA = 1.0*numIgEIgAClones/numIgEClones
	frac_IgEIgGIgA = 1.0*numIgEIgGIgAClones/numIgEClones
	
	printout = dataset + "," 
	printout += str(frac_IgE) + "," + str(frac_IgEIgG) + ","
	printout += str(frac_IgEIgA) + "," + str(frac_IgEIgGIgA)
	meanf.write(printout+"\n")
	
#printout = "Average,"
#printout += str(1.0*sum(mean_eff_IgGs)/len(mean_eff_IgGs)) + ","
#printout += str(1.0*sum(mean_eff_IgAs)/len(mean_eff_IgAs)) + ","
#printout += str(1.0*sum(mean_IgGs)/len(mean_IgGs)) + ","
#printout += str(1.0*sum(mean_IgAs)/len(mean_IgAs))
#meanf.write(printout+"\n")

outf.close()
meanf.close()
