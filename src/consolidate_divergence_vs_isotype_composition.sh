#!/bin/bash
# Consolidate output files from divergence vs isotype composition analysis
COUNT=0
for DATASET in {01107PB,01207PB,SFAPB}
do
    # Create empty file for consolidated results for DATASET
    CONSOLIDATED_FILE=../results/divergence_vs_composition_$DATASET.csv
    touch $CONSOLIDATED_FILE
    
    # For each result file:
    for RESULT_FILE in ../results/${DATASET}_clone_*_divergence_vs_composition_temp.csv
    do
        if [ $COUNT -eq 0 ]; then
            # For the first file, add the first line as the header
            cat $RESULT_FILE >> $CONSOLIDATED_FILE
            #echo -e "\n" >> $CONSOLIDATED_FILE
            COUNT=1
        else
            tail -1 $RESULT_FILE >> $CONSOLIDATED_FILE
            #echo -e "\n" >> $CONSOLIDATED_FILE
        fi
        
        # Remove temporary file
        rm $RESULT_FILE
            
    done

done