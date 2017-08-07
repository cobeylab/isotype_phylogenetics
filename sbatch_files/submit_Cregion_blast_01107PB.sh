#!/bin/bash

# Get all clone files identifiers, concatenate into single variable with ',' for sbatch array:
CLONE_IDS=$(ls ../results/clones/01107PB_clone_* | grep -v 'macse' | grep -o 'clone_[0-9]*'|tr -d [a-z_] | tr '\n' ',')

# Remove trailing comma
CLONE_IDS=${CLONE_IDS::-1}

# Run sbatch command passing clone ids to --array option
sbatch --array=$CLONE_IDS Cregion_blast_01107PB.sbatch
