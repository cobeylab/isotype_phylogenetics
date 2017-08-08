#!/bin/bash

# Get all clone identifiers
ls ../results/clones/01207PB_clone_* | grep -v 'macse' | grep -o 'clone_[0-9]*'|tr -d [a-z_] > clone_ids_01207PB_blast.tmp

# Split clone ids in files with 500 ids (so that sbatch commands do not become too long)
split -l 500 clone_ids_01207PB_blast.tmp 'clone_ids_01207PB_blast' -a 1
rm clone_ids_01207PB_blast.tmp

# Loop over aggregated clone id files, submitting array jobs
for f in clone_ids_01207PB_blast*
do
    #Get clone identifiers, concatenate into single variable with ',' for sbatch array:
    CLONE_IDS=$(cat $f | tr '\n' ',')
    
    # Remove trailing comma
    CLONE_IDS=${CLONE_IDS::-1}

    # Run sbatch command passing clone ids to --array option
    sbatch --array=$CLONE_IDS Cregion_blast_01207PB.sbatch
done

# Remove temporary files
rm clone_ids_01207PB_blast*