#!/bin/bash

# Get all clone identifiers (basically list fasta files excluding macse and temp files)
ls ../results/clones/01107PB_clone_*fasta | grep -v 'macse' | grep -v 'TEMP'| grep -o 'clone_[0-9]*'|tr -d [a-z_] > clone_ids_01107PB_blast.tmp

# Split clone ids in files with 500 ids (so that sbatch commands do not become too long)
split -l 500 clone_ids_01107PB_blast.tmp 'clone_ids_01107PB_blast' -a 1
rm clone_ids_01107PB_blast.tmp

# Loop over aggregated clone id files, submitting array jobs
for f in clone_ids_01107PB_blast*
do
    #Get clone identifiers, concatenate into single variable with ',' for sbatch array:
    CLONE_IDS=$(cat $f | tr '\n' ',')
    
    # Remove trailing comma
    CLONE_IDS=${CLONE_IDS::-1}

    # Run sbatch command passing clone ids to --array option
    sbatch --array=$CLONE_IDS Cregion_blast_01107PB.sbatch
    
    # Wait 120 seconds until new submission
    #sleep 120
done

# Remove temporary files
rm clone_ids_01107PB_blast*