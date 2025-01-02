#!/bin/bash 

# compute_ld_matrix_single_block.sh
# Computes LD matrix for a single LD block using the reference data,
# processes it to reduce file size for SuSiE.

# Load necessary modules (adjust as needed for your environment)
#module load plink/2.0
#module load plink/1.9
#module load bedtools
#module load R  # Load R module for processing LD matrix

# Set paths to directories and files
PGEN_DIR="mounted-data-readonly/merged_arrays"          # Update this to your actual path
LD_BLOCKS_FILE="ld_blocks/ld_blocks_with_ids.bed"         # LD blocks file with block IDs in the 4th column
OUTPUT_DIR="ld_matrices"
mkdir -p "$OUTPUT_DIR"

# Get the block ID from the command line arguments
BLOCK_ID=$1

if [ -z "$BLOCK_ID" ]; then
    echo "No block ID provided. Usage: ./compute_ld_matrix_single_block.sh <BLOCK_ID>"
    exit 1
fi

echo "Computing LD matrix for block $BLOCK_ID"

# Create a temporary directory for this block
TEMP_DIR="temp_ld_${BLOCK_ID}"
mkdir -p "$TEMP_DIR"

# Get the chromosome number for this block
CHR=$(awk -v block="$BLOCK_ID" '$4 == block {print $1}' "$LD_BLOCKS_FILE" | head -1 | sed 's/chr//')

if [ -z "$CHR" ]; then
    echo "Chromosome for block $BLOCK_ID not found. Skipping."
    rm -rf "$TEMP_DIR"
    exit 1
fi

PVAR_FILE="${PGEN_DIR}/gsa_pmra_imputed_chr${CHR}_maf01_rsq8_geno01_ipn_id_excluded_array_assoc.pvar"

if [ ! -f "$PVAR_FILE" ]; then
    echo "PVAR file for chromosome $CHR not found at $PVAR_FILE. Skipping block $BLOCK_ID."
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Extract SNP positions and IDs from the reference data
awk -v chr="$CHR" 'BEGIN{OFS="\t"} !/^#/ {print "chr"chr, $2-1, $2, $3}' "$PVAR_FILE" > "${TEMP_DIR}/chr${CHR}_snps.bed"

# Assign SNPs to LD blocks
bedtools intersect -a "${TEMP_DIR}/chr${CHR}_snps.bed" -b "$LD_BLOCKS_FILE" -wa -wb > "${TEMP_DIR}/chr${CHR}_snps_in_ldblocks.bed"

# Extract SNPs for the block
awk -v block="$BLOCK_ID" '$8 == block {print $4}' "${TEMP_DIR}/chr${CHR}_snps_in_ldblocks.bed" > "${TEMP_DIR}/snps_block_${BLOCK_ID}.txt"

# Check if SNPs exist
if [ ! -s "${TEMP_DIR}/snps_block_${BLOCK_ID}.txt" ]; then
    echo "No SNPs found for block $BLOCK_ID. Skipping."
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Convert PGEN to BED using PLINK 2.0
plink2 --pfile "${PGEN_DIR}/gsa_pmra_imputed_chr${CHR}_maf01_rsq8_geno01_ipn_id_excluded_array_assoc" \
       --extract "${TEMP_DIR}/snps_block_${BLOCK_ID}.txt" \
       --make-bed \
       --out "${TEMP_DIR}/ref_block_${BLOCK_ID}" \
       --threads 1

# Check if extraction was successful
if [ ! -f "${TEMP_DIR}/ref_block_${BLOCK_ID}.bed" ]; then
    echo "Failed to extract genotypes for block $BLOCK_ID. Skipping."
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Compute LD matrix using PLINK 1.9
plink --bfile "${TEMP_DIR}/ref_block_${BLOCK_ID}" \
      --r square \
      --out "${TEMP_DIR}/ld_matrix_block_${BLOCK_ID}" \
      --threads 1

# Process LD matrix to reduce file size
# Run R script to read LD matrix, convert to single precision, save as RDS, and compress
Rscript process_ld_matrix.R "${TEMP_DIR}/ld_matrix_block_${BLOCK_ID}.ld" "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}_single.rds.gz"

# Save SNP IDs (from the .bim file)
cp "${TEMP_DIR}/ref_block_${BLOCK_ID}.bim" "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.bim"

# Remove the temporary LD matrix text file
rm "${TEMP_DIR}/ld_matrix_block_${BLOCK_ID}.ld"

# Clean up temporary files
rm -rf "$TEMP_DIR"

echo "LD matrix for block $BLOCK_ID computed and processed successfully."


