#!/bin/bash
#SBATCH --job-name=run_partis_SFAPB
#SBATCH --output=run_partis_SFAPB.out
#SBATCH --error=run_partis_SFAPB.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

module load python/2.7.12

/project/cobey/partis/partis-0.13.0/bin/partis partition --n-procs 16 --infname /project/cobey/mvieira/isotype_phylogenetics/data/SFAPB_processed.fasta --outfname /project/cobey/mvieira/isotype_phylogenetics/results/SFAPB_partition.csv
