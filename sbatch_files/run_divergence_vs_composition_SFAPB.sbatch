#!/bin/bash
#SBATCH --job-name=run_divergence_vs_composition_SFAPB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o run_divergence_vs_composition_SFAPB_%A_%a.out        
#SBATCH -e run_divergence_vs_composition_SFAPB_%A_%a.err      
#SBATCH --time=200:00:00
#SBATCH --mem-per-cpu=8000

module load python

TREE_FILE=../results/trees/SFAPB_clone_${SLURM_ARRAY_TASK_ID}_annotated_tree.nex

python ../src/divergence_vs_isotype_composition.py $TREE_FILE
