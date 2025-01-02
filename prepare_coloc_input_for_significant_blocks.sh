#!/bin/bash 

# prepare_coloc_input-for_significant_blocks.sh
# Script to prepare coloc input files for significant LD blocks

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a directory and a file name for a GWAS summary statistic."
    echo "Usage: $0 GWAS_SUM_DIR GWAS_SUM_file_name"
    exit 1
fi

# Input and output files
SUMSTATS_DIR=$1
GWAS_FILE_NAME=$2            # Path to your gzipped GWAS summary statistics file
EQTL_DIR=$3
batch=$4

# Input Arguments
GWAS_A_SUMSTATS=${SUMSTATS_DIR}/${GWAS_FILE_NAME} # Path to GWAS A summary statistics
GWAS_B_SUMSTATS_DIR=${EQTL_DIR} # Directory containing GWAS B summary statistics (e.g., gwasB_1.txt, gwasB_2.txt, ..., gwasB_N.txt)
TMP_OUTPUT_DIR="tmp_output"
SIGNIFICANT_BLOCKS="${TMP_OUTPUT_DIR}/significant_ld_blocks.txt"
LD_MATRICES_DIR="filesystems/ld_matrices"
OUTPUT_DIR="coloc_inputs/${batch}" # Directory where coloc input files will be stored

mkdir -p $OUTPUT_DIR

# Delimiter: space-delimited
DELIMITER=" "      # Space delimiter

# Column numbers (adjust if necessary)
CHROM_COLUMN=1    # CHROM
GENPOS_COLUMN=2   # GENPOS
SNPID_COLUMN=3    # ID

# Loop through each significant block and prepare coloc input files
while read -r BLOCK_LINE; do
    # Extract block information
    CHR=$(echo "$BLOCK_LINE" | awk '{print $1}' | sed 's/chr//')
    START=$(echo "$BLOCK_LINE" | awk '{print $2}')
    END=$(echo "$BLOCK_LINE" | awk '{print $3}')
    BLOCK_ID=$(echo "$BLOCK_LINE" | awk '{print $4}')
    echo "Preparing coloc input for block $BLOCK_ID on chromosome $CHR"
  
  # Extract summary statistics for SNPs in the block
    gunzip -c "$GWAS_A_SUMSTATS" | awk -v chr_col=$CHROM_COLUMN -v pos_col=$GENPOS_COLUMN -v snp_col=$SNPID_COLUMN \
        -v chr="$CHR" -v start="$START" -v end="$END" -v FS="$DELIMITER" -v OFS="\t" ' NR==1 || ($chr_col == chr && $pos_col >= start && $pos_col <= end)' \
         > "${OUTPUT_DIR}/gwasA_block_${BLOCK_ID}.txt"

  # Extract relevant summary statistics for the LD block from each GWAS B
  for gwas_b_file in $GWAS_B_SUMSTATS_DIR/ALL-GWAS_*.regenie; do
    gwas_b_id=$(basename $gwas_b_file | sed 's/ALL-GWAS_\(ENSG[0-9]*\)\.regenie/\1/')

    # Define the expected output LD matrix file
    gwasB_BLOCK_FILE="$OUTPUT_DIR/ALL-GWAS_${gwas_b_id}_block_${BLOCK_ID}.txt"
    
    if [ -f "$gwasB_BLOCK_FILE" ]; then
        echo "block $BLOCK_ID already exists at $gwasB_BLOCK_FILE. Skipping creation."
        continue  # Skip to the next LD block
    fi

    {
    # Print the header from the second row
       #gunzip -c "$gwas_b_file" | 
       awk -v FS="$DELIMITER" 'NR == 2' "$gwas_b_file"

       #gunzip -c "$gwas_b_file" | 
       awk -v chr_col=$CHROM_COLUMN -v pos_col=$GENPOS_COLUMN -v snp_col=$SNPID_COLUMN \
        -v chr="$CHR" -v start="$START" -v end="$END" -v FS="$DELIMITER" -v OFS="\t" ' NR>2 && ($chr_col == chr && $pos_col >= start && $pos_col <= end)' "$gwas_b_file"
     } > "$OUTPUT_DIR/ALL-GWAS_${gwas_b_id}_block_${BLOCK_ID}.txt"
  done

  # Copy the corresponding LD matrix for this block
  # cp "$LD_MATRICES_DIR/ld_matrix_block_${block_id}.txt" "$OUTPUT_DIR/ld_matrix_block_${block_id}.txt"

  echo "Prepared coloc input files for Block ID: $block_id"
done < $SIGNIFICANT_BLOCKS

echo "All coloc input files prepared."

