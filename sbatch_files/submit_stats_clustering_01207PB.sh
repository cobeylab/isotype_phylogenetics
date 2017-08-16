ls ../results/trees/01207PB_clone_*annotated_tree.nex | grep -o 'clone_[0-9]*' |tr -d [a-z_] > clone_ids_01207PB_ct.tmp

split -l 500 clone_ids_01207PB_ct.tmp 'clone_ids_01207PB_ct' -a 1
rm clone_ids_01207PB_ct.tmp

for f in clone_ids_01207PB_ct*
do 
	CLONE_IDS=$(cat $f | tr '\n', ',')

	CLONE_IDS=${CLONE_IDS::-1}

	sbatch --array=$CLONE_IDS stats_clustering_01207PB.sbatch

	sleep 120
done

rm clone_ids_01207PB_ct*
