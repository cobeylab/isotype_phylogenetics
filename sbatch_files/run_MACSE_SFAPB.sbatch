#!/bin/bash
#SBATCH --job-name=align_SFAPB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o align_SFAPB_%A_%a.out        
#SBATCH -e align_SFAPB_%A_%a.err        
#SBATCH --time=200:00:00

module load java

MACSE_PATH=/project/cobey/MACSE/macse_git_agro_v1.02/jar_file/macse.jar

INPUT_FILE=SFAPB_clone_${SLURM_ARRAY_TASK_ID}.fasta

# MACSE overrides those output file names for some reason...
OUTPUT_FILE_NT=SFAPB_clone_${SLURM_ARRAY_TASK_ID}_alignment_NT.fasta
OUTPUT_FILE_AA=SFAPB_clone_${SLURM_ARRAY_TASK_ID}_alignment_AA.fasta


java -jar -Xmx600m $MACSE_PATH -prog alignSequences -seq ../results/clones/$INPUT_FILE -../results/clones/$OUTPUT_FILE_NT -../results/clones/$OUTPUT_FILE_AA

# Remove amino acid alignment file
rm ../results/clones/SFAPB_clone_${SLURM_ARRAY_TASK_ID}_macse_AA.fasta
    
# Replace frameshift symbols ('!') from MACSE with gaps ('-')
cat ../results/clones/SFAPB_clone_${SLURM_ARRAY_TASK_ID}_macse_NT.fasta | tr '!' '-' > ${OUTPUT_FILE_NT}.tmp
mv ${OUTPUT_FILE_NT}.tmp ../results/clones/SFAPB_clone_${SLURM_ARRAY_TASK_ID}_macse_NT.fasta




