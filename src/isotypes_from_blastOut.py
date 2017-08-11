import csv
import sys
import glob
from Bio import SeqIO
dataDir = '../results/clones/'

def get_grouped_query_ids(bo_file_name):
	#group same query_id and save them into 'grouped_query_ids'
	bo = open(bo_file_name, "r")
	bo_csv = csv.DictReader(bo)

	query_ids = set()
	for row in bo_csv:
		query_ids.add(row['query_id'])
	query_ids = list(query_ids)

	bo = open(bo_file_name, "r")
	bo_csv = csv.DictReader(bo)

	grouped_query_ids = [[] for i in range(len(query_ids))]
	for row in bo_csv:
		for i in range(len(query_ids)):
			if query_ids[i] == row['query_id']:
				grouped_query_ids[i].append(row)
				break
	bo.close()
	return grouped_query_ids

def get_query_ids_from_fasta(bo_file_name):	
	fname = bo_file_name.replace("_blast_output.csv", ".fasta")
	fastaf = open(fname, "r")
	fasta_object = list(SeqIO.parse(fastaf, "fasta"))
	query_ids_fasta = []
	for obj in fasta_object:
		query_ids_fasta.append(obj.id)
	query_ids_fasta = query_ids_fasta[1:] #excluding "Naive"
	return query_ids_fasta
		

#save file names end with *blast_output.csv within the directory to bo_file_names 
bo_file_names = []
for fnames in glob.glob(dataDir+"*blast_output.csv"):
	bo_file_names.append(fnames)
print (len(bo_file_names))

#write isotype results for each clone 
for bo_file_name in bo_file_names:

	#grouping the same query_id for convenience
	grouped_query_ids = get_grouped_query_ids(bo_file_name)
	query_ids_fasta = get_query_ids_from_fasta(bo_file_name)
	
	#output file for inferred isotypes
	isotype_file_name = bo_file_name.replace("blast_output", "isotypes")
	isotype_file = open(isotype_file_name, "w")
	isotype_file.write("id,isotype,evalue,bit_score,q_start,q_end,s_start,s_end\n")
	
	for same_query_id in grouped_query_ids:
		best_scored = [] #save best scores from blast results for the same query_id

		for one_score in same_query_id:
			id = one_score['query_id']
			qstart = int(one_score['query_start'])
			qend = int(one_score['query_end'])
			sstart = int(one_score['subject_start'])
			send = int(one_score['subject_end'])
			eval = float(one_score['evalue'])
			bitscore = float(one_score['bit_score'])

			#exclude reverse-complement match
			if bool(qstart > qend) != bool(sstart > send) : #exclusive or
				continue				
		
			#exclude if evalue is not low enough
			if eval > 0.01:
				continue
						
			#if there is no element in best_scored, put this score as the best
			if len(best_scored) == 0:
				best_scored.append(one_score)
				continue			
			
			#update best_scored
			if eval < float(best_scored[0]['evalue']):
				best_scored = [one_score]
			elif eval == float(best_scored[0]['evalue']) and bitscore > float(best_scored[0]['bit_score']):
				best_scored = [one_score]
			elif eval == float(best_scored[0]['evalue']) and bitscore == float(best_scored[0]['bit_score']):		
				best_scored.append(one_score)


		#Now decide isotype
		if len(best_scored) == 0:
			isotype = 'Undetermined'
			#we can print any result of blast for this id: it is just either reverse-complement or high evalue
			score_to_print = same_query_id[0] 
			
		elif len(best_scored) > 1:
			#When there are more than one best score, check if the "general isotype"s (IgG or IgA) are the same 			
			general_isotype = best_scored[0]['subject_id'][0:3] #get info only upto IgG, IgA
			same_general_isotype = True
			for other_best in best_scored[1:]:
				if general_isotype != other_best['subject_id'][0:3]:
					same_general_isotype = False
					break
				
			if same_general_isotype == True:
				isotype = general_isotype
			else:
				isotype = 'Ambiguous'
				
			score_to_print = best_scored[0]

		elif len(best_scored) == 1:
			isotype = best_scored[0]['subject_id']
			score_to_print = best_scored[0]

		
		#printout isotype to output file
		printout = ''
		printout += score_to_print['query_id'] +","
		printout += isotype + ","
		printout += score_to_print['evalue'] + ","
		printout += score_to_print['bit_score'] + ","
		printout += score_to_print['query_start'] + ","
		printout += score_to_print['query_end'] + ","
		printout += score_to_print['subject_start'] + ","
		printout += score_to_print['subject_end'] + "\n"

		isotype_file.write(printout)
		
		#query_ids_fasta minus query_ids_blast_output
		#left over ids in query_ids_fasta means they were not blasted because their C-regions were too short 
		query_ids_fasta.pop(query_ids_fasta.index(score_to_print['query_id']))

	#Write results as "Undetermined" for queries that could not be blasted due to short C regions 		
	for short_seq_ids in query_ids_fasta:
		isotype_file.write(short_seq_ids + ",Undetermined,,,,,,\n")

	isotype_file.close()

