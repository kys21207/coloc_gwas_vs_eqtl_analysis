#!/bin/bash

# Ensure both input files are provided
if [ $# -eq 0  ]; then
    echo "Usage: $0 significant_ld_block.txt eqtl.all.tsv.gz eqtl"
    exit 1
fi

SIGNIFICANT_LD_BLOCK=$1
EQTL_FILE=$2
eqtl_name_prefix=$3  
OUTPUT_FILE="tmp_input/${eqtl_name_prefix}_subset.tsv.gz"

# Create a temporary directory for intermediate files
#TMP_DIR=$(mktemp -d)
TMP_DIR=tmp
mkdir -p $TMP_DIR

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Convert significant_ld_block.txt to BED format and sort it
# To change chr1 -> 1
awk 'BEGIN {OFS="\t"} {print substr($1, 4), $2, $3, $0}' $SIGNIFICANT_LD_BLOCK | sort -k1,1 -k2,2n > $TMP_DIR/significant_ld_block.bed

# Extract the header from the eQTL file
#zcat $EQTL_FILE | head -n 1 > $TMP_DIR/header

# Add the header to the temporary directory
echo -e "ensgid\tchr\tpos\tref\talt\tmaf\tbeta\tse" > $TMP_DIR/header

# Convert the eQTL file to BED format
# Only for the eqtl data that has the following schema;  
# molecular_trait_id      chromosome      position        ref     alt     variant ma_samples      maf     pvalue  beta    se      type    ac      an      r2      molecular_trait_object_id       gene_id median_tpm       rsid
zcat $EQTL_FILE | tail -n +2 | awk 'BEGIN {OFS="\t"} {print $2, $3, $3, $1, $2, $3, $4, $5, $8, $10, $11}' | gzip > $TMP_DIR/eqtl.bed.gz
#gzip $TMP_DIR/eqtl.bed

# Use bedtools to intersect the BED files and extract relevant lines
bedtools intersect -a $TMP_DIR/eqtl.bed.gz -b $TMP_DIR/significant_ld_block.bed -wa | cut -f4- > $TMP_DIR/filtered_eqtl

# Combine the header and the filtered eQTL data
cat $TMP_DIR/header $TMP_DIR/filtered_eqtl | gzip > $OUTPUT_FILE

# Clean up temporary files
#rm -r $TMP_DIR

echo "Subset data has been saved to $OUTPUT_FILE"
