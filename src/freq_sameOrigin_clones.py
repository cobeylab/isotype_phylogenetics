import csv
import glob
dataDir = '../results/clones/'

outf = open(dataDir+"frequency_of_sameOrigin_clones.txt", "w")
outf.write("dataset,IgGIgA_together,IgG_noIgA,IgA_noIgG\n")

datasets = ["01107", "01207", "SFA"]
for dataset in datasets:
	total_clones = 0
	sameOrigin_clones = 0
	IgG_noIgA_clones = 0
	IgA_noIgG_clones = 0
	
	for clone_isotypes in glob.glob(dataDir+dataset+"*isotypes.csv"):
		isotype_file = open(clone_isotypes, "r")
		isotype_csv = csv.DictReader(isotype_file)

		isotypes = []
		for sequence in isotype_csv:
			isotypes.append(sequence['isotype'][0:3])
		
		total_clones += 1
				
		if 'IgG' in isotypes and 'IgA' in isotypes:
			sameOrigin_clones += 1
		
		elif 'IgG' in isotypes and 'IgA' not in isotypes:
			IgG_noIgA_clones += 1
		
		elif 'IgA' in isotypes and 'IgG' not in isotypes:
			IgA_noIgG_clones += 1
	
		isotype_file.close()
		
	sameOrigin_freq = 1.0*sameOrigin_clones/total_clones
	IgG_noIgA_freq = 1.0*IgG_noIgA_clones/total_clones
	IgA_noIgG_freq = 1.0*IgA_noIgG_clones/total_clones
	
	printout = dataset + ","
	printout += str(sameOrigin_freq) + ","
	printout += str(IgG_noIgA_freq) + ","
	printout += str(IgA_noIgG_freq) +"\n"
	outf.write(printout)
	
outf.close()
