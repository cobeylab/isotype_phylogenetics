#Get statistics for testing if the frequency of clones having both IgG and IgA is different from random expectation

import csv
import pandas
import numpy

dataDir = "../results/clones/"

obsf = open(dataDir+"frequency_of_sameOrigin_clones.csv", "rU")
obs_csv = csv.DictReader(obsf)
obs_list = []
for row in obs_csv:
	obs_list.append(row)

outf = open(dataDir+"stats_test_freq_sameOrigin_random.csv", "w")
outf.write("dataset,IgGIgA_together_obs,mean_IgGIgA_together_rdm,sd_IgGIgA_together_rdm,CI_min_IgGIgA_together,CI_max_IgGIgA_together,")
outf.write("IgG_noIgA_obs,mean_IgG_noIgA_rdm,sd_IgG_noIgA_rdm,CI_min_IgG_noIgA,CI_max_IgG_noIgA")
outf.write("IgA_noIgG_obs,mean_IgA_noIgG_rdm,sd_IgA_noIgG_rdm,CI_min_IgA_noIgG,CI_max_IgA_noIgG\n")

datasets = ["01107", "01207", "SFA"]
colnames = ['dataset', 'IgGIgA_together', 'IgG_noIgA', 'IgA_noIgG']
for dataset in datasets:

	data = pandas.read_csv(dataDir+'frequency_sameOrigin_random_'+dataset+".csv", names = colnames)
	#[1:] is to remove the header
	dset = data.dataset.tolist()[1:]
	IgGIgA = map(float, data.IgGIgA_together.tolist()[1:])
	IgG_noIgA = map(float, data.IgG_noIgA.tolist()[1:])
	IgA_noIgG = map(float, data.IgA_noIgG.tolist()[1:])
	
	mean_IgGIgA = numpy.mean(IgGIgA)
	sd_IgGIgA = numpy.std(IgGIgA, ddof = 1)
	mean_IgG_noIgA = numpy.mean(IgG_noIgA)
	sd_IgG_noIgA = numpy.std(IgG_noIgA, ddof = 1)
	mean_IgA_noIgG = numpy.mean(IgA_noIgG)
	sd_IgA_noIgG = numpy.std(IgA_noIgG, ddof = 1)
	
	for obs in obs_list:
		if obs['dataset'] == dataset:
			observed = obs
			break
		
	oneline = dataset+","
	obs_IgGIgA = round(float(observed['IgGIgA_together']), 4)
	mean_IgGIgA = round(mean_IgGIgA, 4)
	sd_IgGIgA = round(sd_IgGIgA, 4)
	obs_IgG_noIgA = round(float(observed['IgG_noIgA']), 4)
	mean_IgG_noIgA = round(mean_IgG_noIgA, 4)
	sd_IgG_noIgA = round(sd_IgG_noIgA, 4)
	obs_IgA_noIgG = round(float(observed['IgA_noIgG']), 4)
	mean_IgA_noIgG = round(mean_IgA_noIgG, 4)
	sd_IgA_noIgG = round(sd_IgA_noIgG, 4)

	IgGIgA = sorted(IgGIgA)
	IgG_noIgA = sorted(IgG_noIgA)
	IgA_noIgG = sorted(IgA_noIgG)
	
	alpha = 0.05
	idx = int((alpha/2)*len(IgGIgA))
	IgGIgA_CI = ( round((IgGIgA[idx-1]), 4), round((IgGIgA[-1*idx]), 4) )	
	IgG_noIgA_CI = ( round((IgG_noIgA[idx-1]), 4), round((IgG_noIgA[-1*idx]), 4) )	
	IgA_noIgG_CI = ( round((IgA_noIgG[idx-1]), 4), round((IgA_noIgG[-1*idx]), 4) )	

	oneline += str(obs_IgGIgA)+","+str(mean_IgGIgA)+","+str(sd_IgGIgA)+","+str(IgGIgA_CI[0])+","+str(IgGIgA_CI[1])+","
	oneline += str(obs_IgG_noIgA)+","+str(mean_IgG_noIgA)+","+str(sd_IgG_noIgA)+","+str(IgG_noIgA_CI[0])+","+str(IgG_noIgA_CI[1])+","
	oneline += str(obs_IgA_noIgG)+","+str(mean_IgA_noIgG)+","+str(sd_IgA_noIgG)+","+str(IgA_noIgG_CI[0])+","+str(IgA_noIgG_CI[1])
	
	outf.write(oneline+"\n")

outf.close()
