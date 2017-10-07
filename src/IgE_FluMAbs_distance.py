#IgE_FluMAbs_distance.py: 
#Calculate distance between IgE sequences and FluMbs within same clones and pick the closest IgE to each MAb.

#Codon table dictionary for translation
bases = ['T','C','A','G']
codons = [a+b+c for a in bases for b in bases for c in bases]
aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons,aas))


#Define functions
def translate(nt_seq):
	aa_seq = ''
	for i in range(0, len(nt_seq), 3):
		codon = nt_seq[i:i+3]
		try:
			aa = codon_table[codon]
		except KeyError:
			aa = "?"
		aa_seq += aa
	return aa_seq

def get_best_IgE(distances):
	#Among closest IgEs, get the sequence with longest compared length

	idx_aa_distance = 1
	idx_comp_length = 3
	
	distances.sort(key=lambda x:x[idx_aa_distance])
	best = distances[0]
	for distance in distances:
		if distance[idx_aa_distance] == best[idx_aa_distance]:
			if distance[idx_comp_length] > best[idx_comp_length] :
				best = distance
	return best

def get_most_similar_IgE(MAbseq, IgEseqs):
	#get closest IgE sequence to the MAbseq in terms of amino acid difference

	distances = []
	
	MAb_nt = MAbseq.split("\n")[1]
	MAb_aa = translate(MAb_nt)
	
	for IgEseq in IgEseqs:
		distance = 0
		distance_aa = 0
		lenComp = 0
		
		IgE_nt = IgEseq.split("\n")[1]
		IgE_aa = translate(IgE_nt)
		IgE_name = IgEseq.split("\n")[0].split(">")[1]

		#skip an IgE sequence if its aligned length is different from MAb's length
		if len(MAb_nt) != len(IgE_nt):
			continue

		#nucleotide distance
		for s in range(len(MAb_nt)):
			#Compare nt bases only at the site where both sequences have nt base
			if MAb_nt[s] == 'N' or IgE_nt[s] == 'N':
				continue
			else:
				if MAb_nt[s] != IgE_nt[s]:
					distance += 1
		#amino acid distance
		for s in range(len(MAb_aa)):
			#Compare AA only at the site where both sequences have AA
			if MAb_aa[s] == '?' or IgE_aa[s] == '?':
				continue
			else:
				lenComp += 1
				if MAb_aa[s] != IgE_aa[s]:
					distance_aa += 1
		
		distances.append((IgE_name, distance_aa, distance, lenComp, IgE_aa, IgE_nt))
		
	best = get_best_IgE(distances)
	
	return best

def read_original_fasta(originalf):
	originals = []
	for line in originalf:
		if line.find(">") >= 0:
			seq = line
		else:
			seq += line
			originals.append(seq)
	return originals
	
def put_seq_in_fasta(IgEname, distance_fastaf, originals):
	for seq in originals:
		if IgEname in seq:
			distance_fastaf.write(seq)
			return IgEname
			break
	return 0

#file names
fastaf_path = "../results/clones/"

originalf_name = "../data/SFAPB.fasta"
isotypef_name = "../results/SFAPB_FluMAbs_clone_isotypes.csv"
distancef_name = "../results/SFAPB_FluMAbs_closeIgE_byAA.csv"
pickedf_name = "../results/SFAPB_FluMAbs_IgE_picked.csv"
distance_fastaf_name = "../results/SFAPB_FluMAbs_closeIgE_byAA.fasta"
picked_fastaf_name = "../results/SFAPB_FluMAbs_IgE_picked.fasta"

#open files and write head
originalf = open(originalf_name, "rU")
isotypef = open(isotypef_name, "r")
distancef = open(distancef_name, "w")
pickedf = open(pickedf_name, "w")
distance_fastaf = open(distance_fastaf_name, "w")
picked_fastaf = open(picked_fastaf_name, "w")

distancef.write("clone_name,mAb_sequence,closest_IgE,distance_aa,distance_nt,len_aa_compared,num_mAbs_in_clone,num_IgE_in_clone\n")
pickedf.write("clone_name,mAb_sequence,closest_IgE,distance_aa,distance_nt,len_aa_compared,num_mAbs_in_clone,num_IgE_in_clone\n")

# read original fasta file
originals = read_original_fasta(originalf)

#lists
distance_list = []

clone_list = []
mAb_list = []
seqnames_list = []
isotypes_list = []

#put clones, mAbs, seqs, and isotypes in lists, indexing matched to clone_list
for line in isotypef:
	each = line.split("\n")[0].split(",")
	mAb = each[0]
	clone = each[1]
	seqnames = each[2]
	isotypes = each[3]
	if clone in clone_list:
		idx = clone_list.index(clone)
		mAb_list[idx].append(mAb)
	else:
		clone_list.append(clone)
		mAb_list.append([mAb])
		seqnames_list.append(seqnames.split("|"))
		isotypes_list.append(isotypes.split(" "))

#for each clone, get IgE seuqnece names
IgEseqnames_list = [[] for i in range(len(isotypes_list))]

for i in range(len(isotypes_list)):
	for j in range(len(isotypes_list[i])):
		if isotypes_list[i][j] == "IgE":
			IgEseqnames_list[i].append(seqnames_list[i][j])
	
#for each clone, get MAbs and IgE seqs and calculate distance
for idx in range(len(clone_list)):

	#If clone is not assinged to mAb or there is no IgE in the clone, continue.
	if clone_list[idx] == "Not assigned":
		continue
	if len(IgEseqnames_list[idx]) == 0:
		continue

	clonef_name = fastaf_path + clone_list[idx]+"_withN.fasta"
	clonef = open(clonef_name, "rU")

	#make lists of MAbseqs and IgEseqs for this clone
	IgEseqs = []
	MAbseqs = []
	for line in clonef:
		if len(line) < 2:
			continue
		if line.find(">") >= 0:
			IgEseq = ''
			MAbseq = ''
			seqName = line.split("\n")[0].replace(">", "")
			if seqName in IgEseqnames_list[idx]:
				IgEseq = line
			elif seqName in mAb_list[idx]:
				MAbseq = line
				
		else:
			if IgEseq != '':
				nt_seq = line.split("\n")[0]
				IgEseq += nt_seq
				IgEseqs.append(IgEseq)
			elif MAbseq != '':
				nt_seq = line.split("\n")[0]
				MAbseq += nt_seq
				MAbseqs.append(MAbseq)

	clone_name = clonef_name.split("clones/")[1].split("_withN")[0]

	if len(IgEseqs) == 0:
		continue

	#For each MAbseqs, get the most similar IgEseq
	for MAbseq in MAbseqs:
		MAbname = MAbseq.split("\n")[0].split(">")[1]
 		
		closest =  get_most_similar_IgE(MAbseq, IgEseqs)
		IgEname = closest[0]
		distance_aa = closest[1]
		distance = closest[2]
		lenComp = closest[3]
		aa = closest[4] 
		nt = closest[5]
		
		oneline = clone_name+","+MAbname+","+IgEname+","
		oneline += str(distance_aa)+","+str(distance)+","+str(lenComp)+","
		oneline += str(len(MAbseqs))+","+str(len(IgEseqs))
		distancef.write(oneline+"\n")
		
		distance_list.append((clone_name,MAbname,IgEname,distance_aa,distance,lenComp,len(MAbseqs),len(IgEseqs), aa, nt))		

		put_seq_in_fasta(IgEname, distance_fastaf, originals)

distancef.close()
isotypef.close()

#indices for distance_list
idx_clone = 0
idx_mAb = 1
idx_IgE = 2
idx_aa = 3
idx_nt = 4
idx_len = 5
idx_num_mAb = 6
idx_num_IgE = 7
idx_aaseq = 8
idx_ntseq = 9

good = []
bad = []

by_distance = []
by_MAb_freq = []
by_IgE_freq = []
by_0_and_len = []

#pick one 0-distance IgE from a clone of largest number of mAb
distance_list.sort(key=lambda x:x[idx_num_mAb], reverse=True)
most_mAbs = []
for distance in distance_list:
	if distance[idx_num_mAb] == distance_list[0][idx_num_mAb]:
		most_mAbs.append(distance)
most_mAbs.sort(key=lambda x:x[idx_aa])
for seq in most_mAbs:
	if seq[idx_aa] == 0:
		by_MAb_freq.append(seq)
		break

picked = [str(i) for i in by_MAb_freq[0][:-2]]
pickedf.write(",".join(picked) + "\n")
put_seq_in_fasta(picked[idx_IgE], picked_fastaf, originals)


#pick one 0-distance IgE from clone having largest number of IgE 
distance_list.sort(key=lambda x:x[idx_num_IgE], reverse=True)
most_IgEs = []
for distance in distance_list:
	if distance[idx_num_IgE] == distance_list[0][idx_num_IgE]:
		most_IgEs.append(distance)
for seq in most_IgEs:
	if seq[idx_aa] == 0:		
		by_IgE_freq.append(seq)
picked = [str(i) for i in by_IgE_freq[0][:-2]]
pickedf.write(",".join(picked) + "\n")
put_seq_in_fasta(picked[idx_IgE], picked_fastaf, originals)

#pick IgEs with 1-3 distance. When are are more than 1 such IgE from one clone, choose the best one. 
for distance in distance_list:
	if distance[idx_aa] > 0 and distance[idx_aa] < 4 and distance[idx_len] > 100:
		if len(good) > 0 and distance[idx_clone] == good[-1][idx_clone]:
			if distance[idx_aa] < good[-1][idx_aa]:
				good.pop()
				good.append(distance)
			else:
				continue
		else:
			good.append(distance)
	else:
		bad.append(distance)

for g in good:
	picked = [str(i) for i in g[:-2]]
	pickedf.write( ','.join(picked) + "\n")
	put_seq_in_fasta(picked[idx_IgE], picked_fastaf, originals)
	
