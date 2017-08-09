import csv
import sys
from heapq import nlargest
csv.field_size_limit(sys.maxsize)

dataDir = '../results/'
datasets = ["01107PB", "01207PB", "SFAPB"]

clone_size_allseqs = []
clone_size_prdseqs = []
num_unannotated = 0

for data in datasets:
	partis_file = open(dataDir+data+"_partition-cluster-annotations.csv", "r")
	partis_csv = csv.DictReader(partis_file)
	
	for row in partis_csv:
		in_frame = row['in_frames'].split(":")
		has_stop_codons = row['stops'].split(":")
	
		#excluding sequences for which naive sequences could not be found
		if in_frame == ['']:
			num_unannotated += 1
			continue

		#save the number of sequences within the clone regardless of productive/unproductive
		clone_size_allseqs.append(len(in_frame))
			
		#save the size of the clone in terms of the number of productive sequences
		n_productive_seqs = [i for i in range(len(in_frame)) if in_frame[i] == 'True' and has_stop_codons[i] == 'False']
		n_productive_seqs = len(n_productive_seqs)
		if n_productive_seqs == 0: #it is not clone if all sequences are unproductive
			continue	
		clone_size_prdseqs.append(n_productive_seqs)

print num_unannotated

sizes_prdseqs = sorted(list(set(clone_size_prdseqs)))
sizes_allseqs = sorted(list(set(clone_size_allseqs)))
count_sizes_prdseqs = []
count_sizes_allseqs = []

for size in sizes_prdseqs:
	count_sizes_prdseqs.append(clone_size_prdseqs.count(size))
for size in sizes_allseqs:
	count_sizes_allseqs.append(clone_size_allseqs.count(size))


outf = open(dataDir+"clone_size_freqdist_productive.csv", "w")
outf.write("clone_size, count, frequency\n")
#exclude clones with size 0
total_num_clones = sum(count_sizes_prdseqs)
for i in range(len(sizes_prdseqs)):
	freq = 1.0*count_sizes_prdseqs[i]/total_num_clones
	line =  str(sizes_prdseqs[i]) + ", " + str(count_sizes_prdseqs[i]) + ", "+ str(freq)
	outf.write(line+"\n")
outf.close()

outf = open(dataDir+"clone_size_freqdist_allseqs.csv", "w")
outf.write("clone_size, count, frequency\n")
#exclude clones with size 0
total_num_clones = sum(count_sizes_allseqs)
for i in range(len(sizes_allseqs)):
	freq = 1.0*count_sizes_allseqs[i]/total_num_clones
	line =  str(sizes_allseqs[i]) + ", " + str(count_sizes_allseqs[i]) + ", "+ str(freq)
	outf.write(line+"\n")
outf.close()

