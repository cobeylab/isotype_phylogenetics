#!/bin/bash

# Get clone identifiers for all clones with raxml best trees
ls ../results/trees/RAxML_bestTree.01107PB_clone_* | grep -o 'clone_[0-9]*'|tr -d [a-z_] > clone_ids_01107PB_divergence.tmp

# Split clone ids in files with 500 ids (so that sbatch commands do not become too long)
split -l 500 clone_ids_01107PB_divergence.tmp 'clone_ids_01107PB_annotation' -a 1
rm clone_ids_01107PB_divergence.tmp

# Loop over aggregated clone id files, submitting array jobs
for f in clone_ids_01107PB_annotation*
do
    #Get clone identifiers, concatenate into single variable with ',' for sbatch array:
    CLONE_IDS=$(cat $f | tr '\n' ',')
    
    # Remove trailing comma
    CLONE_IDS=${CLONE_IDS::-1}

    # Run sbatch command passing clone ids to --array option
    sbatch --array=$CLONE_IDS run_divergence_vs_composition_01107PB.sbatch
done

# Remove temporary files
rm clone_ids_01107PB_annotation*
