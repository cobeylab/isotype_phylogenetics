#!/bin/bash
#SBATCH --job-name=run_partis_01107PB
#SBATCH --output=run_partis_01107PB.out
#SBATCH --error=run_partis_01107PB.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

module load python/2.7.12

/project/cobey/partis/partis-0.13.0/bin/partis partition --n-procs 16 --infname /project/cobey/mvieira/isotype_phylogenetics/data/01107PB_processed.fasta --outfname /project/cobey/mvieira/isotype_phylogenetics/results/01107PB_partition.csv
