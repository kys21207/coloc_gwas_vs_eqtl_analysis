#!/bin/bash 

# identify_significant_ld_blocks.sh
# Identifies LD blocks that contain significant SNPs.

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a directory and a file name for a GWAS summary statistic."
    echo "Usage: $0 GWAS_SUM_DIR GWAS_SUM_file_name"
    exit 1
fi

# Input and output files
gwas_name_prefix=$1            # Path to your gzipped GWAS summary statistics file

# Input files
TMP_INPUT_DIR="tmp_input"
TMP_OUTPUT_DIR="tmp_output"
FILTERED_SUMSTATS_FILE="${TMP_INPUT_DIR}/${gwas_name_prefix}_subset.tsv.gz"
LD_BLOCKS_FILE="ld_blocks/ld_blocks_with_ids.bed"

mkdir -p "$TMP_OUTPUT_DIR"

# Output file
SIGNIFICANT_LD_BLOCKS_FILE="${TMP_OUTPUT_DIR}/significant_ld_blocks.txt"

# Delimiter: space-delimited
DELIMITER=" "      # Space delimiter

#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N BETA SE CHISQ LOG10P INFO

# Column numbers (adjust if necessary)
CHROM_COLUMN=1    # CHROM is column 1
GENPOS_COLUMN=2   # GENPOS is column 2
SNPID_COLUMN=3    # ID is column 3

# Convert summary statistics to BED format
gunzip -c "$FILTERED_SUMSTATS_FILE" | awk -v chr_col=$CHROM_COLUMN -v pos_col=$GENPOS_COLUMN -v snp_col=$SNPID_COLUMN -v FS="$DELIMITER" -v OFS="\t" \
    'NR>1 {print "chr"$chr_col, $pos_col-1, $pos_col, $snp_col}'  > "${TMP_OUTPUT_DIR}/significant_snps.bed"

bedtools intersect -a "${LD_BLOCKS_FILE}" -b "${TMP_OUTPUT_DIR}/significant_snps.bed" -wa | sort | uniq > "${SIGNIFICANT_LD_BLOCKS_FILE}"

