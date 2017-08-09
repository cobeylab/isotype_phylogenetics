#!/bin/bash
#SBATCH --job-name=run_raxml_01207PB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o run_raxml_01207PB_%A_%a.out        
#SBATCH -e run_raxml_01207PB_%A_%a.err      
#SBATCH --time=200:00:00

ALIGNMENT_FILE=../results/clones/01207PB_clone_${SLURM_ARRAY_TASK_ID}_macse_NT.fasta
OUTPUT_ID=01207PB_clone_${SLURM_ARRAY_TASK_ID}

/project/cobey/RaxML/raxmlHPC-AVX -s $ALIGNMENT_FILE -n $OUTPUT_ID -m GTRGAMMA -o NAIVE -p 57918

# Moving files to results/trees
mv RAxML_*.$OUTPUT_ID ../results/trees/

# Removing reduced fasta file, if any
rm ../results/clones/01207PB_clone_${SLURM_ARRAY_TASK_ID}_macse_NT.fasta.reduced

    