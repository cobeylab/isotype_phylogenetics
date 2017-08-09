#!/bin/bash

# Align only clones with at least MINSEQ sequences:
MINSEQ=10

# Get all clone identifiers
ls ../results/clones/SFAPB_clone_* | grep -v 'macse' | grep -v '.csv'| grep -o 'clone_[0-9]*'|tr -d [a-z_] | head > clone_ids_SFAPB.tmp

# Split clone ids in files with 500 ids (so that sbatch commands do not become too long)
split -l 500 clone_ids_SFAPB.tmp 'clone_ids_SFAPB' -a 1
rm clone_ids_SFAPB.tmp

# Loop over aggregated clone id files, submitting array jobs
for f in clone_ids_SFAPB*
do
    #Get clone identifiers, concatenate into single variable with ',' for sbatch array:
    CLONE_IDS=$(cat $f | tr '\n' ',')
    
    # Remove trailing comma
    CLONE_IDS=${CLONE_IDS::-1}
    
    # Retain only clones with at least MINSEQ sequences
    RETAINED_IDS=''
    
    # For each clone id
    for ID in ${CLONE_IDS//,/ }
    do
      # Number of sequences in clone
      NUMSEQ=$(grep '>' ../results/clones/SFAPB_clone_${ID}.fasta | wc -l)

      # Naive sequence doesn't count towards minimum number of sequences
      NUMSEQ=$(expr $NUMSEQ - $'1')
      
      # Retain if clone has minimum number of sequences
      if [ "$NUMSEQ" -ge "$MINSEQ" ]
        then
        RETAINED_IDS=${RETAINED_IDS},$ID
      
      fi
 
    done

    # Remove leading comma
    RETAINED_IDS=${RETAINED_IDS:1} 
    
    # Run sbatch command passing retained clone ids to --array option
    sbatch --array=$RETAINED_IDS run_MACSE_SFAPB.sbatch       
    
done

# Remove temporary files
rm clone_ids_SFAPB*