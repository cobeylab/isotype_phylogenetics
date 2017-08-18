ls ../results/trees/SFAPB_clone_*annotated_tree.nex | grep -o 'clone_[0-9]*' |tr -d [a-z_] > clone_ids_SFAPB_ct.tmp

split -l 500 clone_ids_SFAPB_ct.tmp 'clone_ids_SFAPB_ct' -a 1
rm clone_ids_SFAPB_ct.tmp

for f in clone_ids_SFAPB_ct*
do 
	CLONE_IDS=$(cat $f | tr '\n', ',')

	CLONE_IDS=${CLONE_IDS::-1}

	sbatch --array=$CLONE_IDS stats_clustering_SFAPB.sbatch

	sleep 120
done

rm clone_ids_SFAPB_ct*
