#!/bin/bash
#SBATCH --job-name=run_partis_01207PB
#SBATCH --output=run_partis_01207PB.out
#SBATCH --error=run_partis_01207PB.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

module load python/2.7.12
module load mafft/7.310

/project/cobey/partis/partis-0.13.0/bin/partis partition --n-procs 16 --infname /project/cobey/mvieira/isotype_phylogenetics/data/01207PB_processed.fasta --outfname /project/cobey/mvieira/isotype_phylogenetics/results/01207PB_partition.csv
