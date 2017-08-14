#!/bin/bash
#SBATCH --job-name=annotate_trees01107PB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o annotate_trees01107PB_%A_%a.out        
#SBATCH -e annotate_trees01107PB_%A_%a.err      
#SBATCH --time=200:00:00

module load python

TREE_FILE=../results/trees/RAxML_bestTree.01107PB_clone_${SLURM_ARRAY_TASK_ID}

python ../src/annotate_tree.py $TREE_FILE
