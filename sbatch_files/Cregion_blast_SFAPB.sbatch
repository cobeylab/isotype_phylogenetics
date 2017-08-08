#!/bin/bash
#SBATCH --job-name=Cregion_blast_SFAPB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o Cregion_blast_SFAPB_%A_%a.out        
#SBATCH -e Cregion_blast_SFAPB_%A_%a.err        
#SBATCH --time=200:00:00

module load blastplus

../src/determine_clone_isotypes.sh ../results/clones/SFAPB_clone_${SLURM_ARRAY_TASK_ID}.fasta


