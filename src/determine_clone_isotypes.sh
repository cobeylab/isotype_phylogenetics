#!/bin/bash
# Takes a clone, finds its sequence identifiers, looks in the processed fasta files for the original sequences (including the constant regions), and blasts them against c_regions.fasta

CLONE_FILE=$1

# Path to output file with blast results
BLAST_OUTPUT=$(echo $CLONE_FILE | grep -o [A-Z]*[0-9]*PB_clone_[0-9]*)_Cregions.csv
BLAST_OUTPUT=../results/clones/$BLAST_OUTPUT

# Path to file with constant region sequences
CREGION_FASTA=/project/cobey/mvieira/isotype_phylogenetics/src/c_regions.fasta

# Find processed fasta file w/ full sequences from dataset identifier
DATASET_ID=$(echo $CLONE_FILE | grep -o [A-Z]*[0-9]*PB)
FULLSEQ_FASTA=../data/${DATASET_ID}_processed.fasta

# Temporary fasta file with original sequences but removed new lines (for grep)
FULLSEQ_TEMP=../results/clones/$(echo $CLONE_FILE | grep -o '[A-Z]*[0-9]*PB_clone_[0-9]*')_TEMP_FULLSEQ.fasta
cat $FULLSEQ_FASTA | tr -d '\n' > $FULLSEQ_TEMP

# Temporary fasta file with query sequences
TEMP_FASTA=$(echo $CLONE_FILE | grep -o [A-Z]*[0-9]*PB_clone_[0-9]*)_TEMP.fasta
TEMP_FASTA=../results/clones/$TEMP_FASTA

touch $TEMP_FASTA

# Get all sequence identifiers, concatenate into comma-separated string
SEQ_IDS=$(cat $CLONE_FILE | grep '>' | grep -v 'NAIVE' | tr '\n' ',')

# Remove trailing comma:
SEQ_IDS=${SEQ_IDS::-1}

# For each sequence identifier in clone fasta file
for ID in ${SEQ_IDS//,/ }
    do
    
    # Find full sequence (including constant region) from original fasta file:
    FULL_SEQ=$(grep -a -o "$ID[A-Z]*" $FULLSEQ_TEMP)
    
    # Remove ID from sequence string:
    FULL_SEQ=${FULL_SEQ/$ID/}
    
    # Write ID and sequence to temporary fasta file:
    #echo $ID
    
    echo $ID >> $TEMP_FASTA
    
    #echo $FULL_SEQ
    echo $FULL_SEQ >> $TEMP_FASTA
    echo '' >> $TEMP_FASTA
    
    done 
       
# Blast $TEMP_FASTA against c_regions.fasta    
blastn -task blastn -query $TEMP_FASTA -out $BLAST_OUTPUT -outfmt 10 -subject $CREGION_FASTA

# Add header to blast output file (from -outfmt 7)
HEADER='query_id,subject_id,percent_identity,alignment_length,mismatches,gap_opens,query_start,query_end,subject_start,subject_end,evalue,bit_score'

echo $HEADER | cat - $BLAST_OUTPUT > ${BLAST_OUTPUT}.tmp && mv ${BLAST_OUTPUT}.tmp $BLAST_OUTPUT

# Remove temporary files
rm $TEMP_FASTA
rm $FULLSEQ_TEMP