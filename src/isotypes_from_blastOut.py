import csv
import sys
import glob

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

#save file names end with *Cregions.csv within the directory
bo_file_names = []
for fnames in glob.glob(dataDir+"*Cregions.csv"):
	bo_file_names.append(fnames)
print (len(bo_file_names))
 
for bo_file_name in bo_file_names:
	print (bo_file_name)
	#grouping the same query_id for convenience
	grouped_query_ids = get_grouped_query_ids(bo_file_name)
	
	#output file for inferred isotypes
	isotype_file_name = bo_file_name.replace("Cregions", "isotypes")
	isotype_file = open(isotype_file_name, "w")
	isotype_file.write("id, isotype, evalue, bit_score, q_start, q_end, s_start, s_end\n")
	
	for same_query_id in grouped_query_ids:
#		print (same_query_id[0]['query_id'])
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
#				print ("  reverse: ", eval, bitscore, qstart, qend, sstart, send)
				continue				
		
			#exclude if evalue is not low enough
			if eval > 0.01:
#				print("  high evalue: ", eval)
				continue
						
			#if there is no element in best_scored, put this score as the best
			if len(best_scored) == 0:
				best_scored.append(one_score)
				print ("  appended to best scored as the first element: ", eval, bitscore)
				continue
				
			
			#update best_scored
			if eval < float(best_scored[0]['evalue']):
				best_scored = [one_score]
#				print ("  replaced due to evalue: ", eval, best_scored[0]['evalue'])
			elif eval == float(best_scored[0]['evalue']) and bitscore > float(best_scored[0]['bit_score']):
				best_scored = [one_score]
#				print ("  replaced due to bitscore: ", bitscore, best_scored[0]['bit_score'])
			elif eval == float(best_scored[0]['evalue']) and bitscore == float(best_scored[0]['bit_score']):		
				best_scored.append(one_score)
#				print ("  appended: ", eval, bitscore, best_scored[0]['evalue'], best_scored[0]['bit_score'])


		#Now decide isotype
		if len(best_scored) == 0:
			isotype = 'Unknown'
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

#		print ("      isotype: ", isotype)
		
		#printout isotype to output file
		printout = ''
		printout += score_to_print['query_id'] +", "
		printout += isotype + ", "
		printout += score_to_print['evalue'] + ", "
		printout += score_to_print['bit_score'] + ", "
		printout += score_to_print['query_start'] + ", "
		printout += score_to_print['query_end'] + ", "
		printout += score_to_print['subject_start'] + ", "
		printout += score_to_print['subject_end'] + "\n"

		isotype_file.write(printout)
		
	isotype_file.close()

