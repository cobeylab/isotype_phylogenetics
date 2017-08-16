#!/bin/bash
#SBATCH --job-name=stats_clustering_01207PB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o stats_clustering_01207PB_%A_%a.out
#SBATCH -e stats_clustering_01207PB_%A_%a.err
#SBATCH --time=200:00:00

module load python

python ../src/stats_for_test_clustering.py ../results/trees/01207PB_clone_${SLURM_ARRAY_TASK_ID}_annotated_tree.nex ../results/test_clustering/stats_clustering_01207PB_clone_${SLURM_ARRAY_TASK_ID}.txt
