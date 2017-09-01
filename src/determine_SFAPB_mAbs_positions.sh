#!/bin/bash
# Takes the sequence identifiers of  SFA flu positive mAbs, identifies the clones each one belongs to, and counts the numbers of the different isotypes in those clones.

FLU_MABS_IDS=$(grep ">" ../data/SFAPB_FluMAbs.fasta | tr -d '>')


# Create summary CSV file
SUMMARY_FILE=../results/SFAPB_FluMAbs_clone_positions.csv
echo "mAb_sequence,clone_name,n_sequences,n_igg,n_iga,n_igm,n_ige,n_igd,n_ambiguous,n_undetermined" > $SUMMARY_FILE

# For each mAb id
for MAB_ID in $FLU_MABS_IDS
do
    echo $MAB_ID
    # Find clone fasta file that contains that mAb id:
    MATCHING_FILE=$(grep -E -l "$MAB_ID" ../results/clones/SFAPB_HT_plus_FluMAbs_clone_*)
    
    # If mAb has been assigned to a clone
    if [ -n "$MATCHING_FILE" ]
    then
    
        # From file name, find clone name:
        CLONE_NAME=$(echo $MATCHING_FILE | grep -o "SFAPB_HT_plus_FluMAbs_clone_[0-9]*")
    
        # Compute clone size
        CLONE_SIZE=$(grep '>' $MATCHING_FILE | wc -l)
    
        # Remove naive sequence from clone size
        CLONE_SIZE=$(($CLONE_SIZE -1 ))
        
        # Get a list of all sequence identifiers in that clone (minus naive and mAb MAB_IDs):
        SEQ_IDS=$(grep '>' $MATCHING_FILE | grep -v '>NAIVE' | grep -v ">007-" | tr -d '>')
    
        # Find the isotypes of all sequences in the clone
        # (| separated list of sequence ids to pass to grep)
        SEQ_IDS_GREP=$(echo $SEQ_IDS | tr ' ' '|')
    
        # If there's at least one sequence in the clone that's not naive nor mAb:
        if [ -n "$SEQ_IDS_GREP" ]
        then
            ISOTYPE_LIST=$(grep -E "$SEQ_IDS_GREP" ../results/clones/SFAPB*isotypes* | cut -d ':' -f 2 | cut -d ',' -f 2)
        else
            ISOTYPE_LIST=''
        fi
        
        #Count isotypes
        N_IGG=$(echo $ISOTYPE_LIST | tr ' ' '\n' | grep 'IgG' | wc -l)
        N_IGA=$(echo $ISOTYPE_LIST | tr ' ' '\n' | grep 'IgA' | wc -l) 
        N_IGM=$(echo $ISOTYPE_LIST | tr ' ' '\n' | grep 'IgM' | wc -l)
        N_IGE=$(echo $ISOTYPE_LIST | tr ' ' '\n' | grep 'IgE' | wc -l)
        N_IGD=$(echo $ISOTYPE_LIST | tr ' ' '\n' | grep 'IgD' | wc -l)
        N_AMBIGUOUS=$(echo $ISOTYPE_LIST | tr ' ' '\n' | grep 'Ambiguous' | wc -l)
        N_UNDETERMINED=$(echo $ISOTYPE_LIST | tr ' ' '\n' | grep 'Undetermined' | wc -l)
    
    # If sequence has not been added to a clone:
    else
        CLONE_NAME='Not assigned'
        CLONE_SIZE='NA'
        N_IGG='NA';N_IGA='NA';N_IGM='NA';N_IGE='NA';N_IGD='NA';N_AMBIGUOUS='NA';N_UNDETERMINED='NA'
    fi
    
    # Add results to summary file:
    echo ${MAB_ID/>/},${CLONE_NAME},${CLONE_SIZE},${N_IGG},${N_IGA},${N_IGM},${N_IGE},${N_IGD},${N_AMBIGUOUS},${N_UNDETERMINED} >> $SUMMARY_FILE
    
done



