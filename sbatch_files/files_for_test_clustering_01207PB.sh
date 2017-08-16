#!/bin/bash
#This file contains summaries of clustering statistics of all clones
SMR_FILE=../results/test_clustering/summary_stats_clustering_01207PB.csv

count=0

for f in $(ls ../results/test_clustering/stats_clustering_01207PB_clone_*txt)
do
	CLONE_ID=$(echo $f | grep -o 'clone_[0-9]*')

	#write the summary line (which is the last line of stats_clustering*.txt file) to the summary file
	#for the first stats_clustering*.txt file, write last 2 lines to the summary file to write a head
	if [ $count = 0 ]; then
		cat $f | tail -n 2 > $SMR_FILE
	else
		cat $f | tail -n 1 >> $SMR_FILE
	fi

	#For each clone, write a file (stats_clustering_*_rdm.csv) that contains 1000 randomized clustering statistics
	RDM_FILE=../results/test_clustering/stats_clustering_01207PB_${CLONE_ID}_rdm.csv
	#except for last 2 lines, write all lines to stats_clustering_*_rdm.csv file
	cat $f | head -n -2 > $RDM_FILE

	count=$((count+1))
done
