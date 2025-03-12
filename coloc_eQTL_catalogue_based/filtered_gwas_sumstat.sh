#!/bin/bash 

# filter_sumstats.sh
# Filters GWAS summary statistics for SNPs with LOG10P > 5 (p-value < 1e-5).
# Handles input files compressed in .gz format and space-delimited.

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a directory and a file name for a GWAS summary statistic."
    echo "Usage: $0 GWAS_SUM_DIR GWAS_SUM_file_name"
    exit 1
fi

# Input and output files
SUMSTATS_FILE=$1            # Path to your gzipped GWAS summary statistics file
gwas_name_prefix=$2 
FILTERED_SUMSTATS_FILE="tmp_input/${gwas_name_prefix}_subset.tsv"

#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N BETA SE CHISQ LOG10P INFO
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N BETA SE CHISQ LOG10P INFO rsID
# Determine if the file is space-delimited by inspecting the first data line
FIRST_LINE=$(gunzip -c "$SUMSTATS_FILE" | grep -v "^CHROM" | head -n 1)

# Count the number of fields in the first data line
NUM_FIELDS=$(echo "$FIRST_LINE" | awk '{print NF}')

echo "Detected $NUM_FIELDS fields in the summary statistics file."

# Set LOG10P_COLUMN based on the number of fields
if [ "$NUM_FIELDS" -eq 14 ]; then
    LOG10P_COLUMN=13
elif [ "$NUM_FIELDS" -eq 12 ]; then
    LOG10P_COLUMN=11
elif [ "$NUM_FIELDS" -eq 13 ]; then
    LOG10P_COLUMN=11
else
    echo "Error: Unexpected number of fields ($NUM_FIELDS). Expected 12 or 14."
    exit 1
fi

echo "LOG10P_COLUMN set to $LOG10P_COLUMN."


# Delimiter: space-delimited
DELIMITER=" "      # Space delimiter
#DELIMITER="\t"      # Tab delimiter

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$FILTERED_SUMSTATS_FILE")"

# Check if the summary statistics file exists
if [ ! -f "$SUMSTATS_FILE" ]; then
    echo "Error: Summary statistics file '$SUMSTATS_FILE' not found."
    exit 1
fi

# Determine if the file is space-delimited by inspecting the first data line
FIRST_LINE=$(gunzip -c "$SUMSTATS_FILE" | grep -v "^CHROM" | head -n 1)

# Count the number of fields in the first data line
NUM_FIELDS=$(echo "$FIRST_LINE" | awk '{print NF}')

echo "Detected $NUM_FIELDS fields in the summary statistics file."

# Verify that LOG10P_COLUMN does not exceed the number of fields
if [ "$LOG10P_COLUMN" -gt "$NUM_FIELDS" ]; then
    echo "Error: LOG10P_COLUMN ($LOG10P_COLUMN) exceeds the number of fields ($NUM_FIELDS)."
    exit 1
fi

# Filter SNPs with LOG10P > 5 and ensure LOG10P is numeric
echo "Filtering SNPs with LOG10P > 5..."

gunzip -c "$SUMSTATS_FILE" | awk -v log10p_col=$LOG10P_COLUMN -v FS="$DELIMITER" -v OFS="\t" '
NR == 1 {print; next} 
{
    # Check if LOG10P is a number and greater than 5
    if ($log10p_col ~ /^-?[0-9]*\.?[0-9]+$/ && $log10p_col > 5) {
        print
    }
}' > "$FILTERED_SUMSTATS_FILE"


# Count the number of SNPs after filtering
TOTAL_SNPS=$(gunzip -c "$SUMSTATS_FILE" | wc -l)
FILTERED_SNPS=$(wc -l < "$FILTERED_SUMSTATS_FILE")
gzip $FILTERED_SUMSTATS_FILE

echo "Total SNPs in original file: $TOTAL_SNPS"
echo "Total SNPs after filtering: $FILTERED_SNPS"

# Optional: Report how many SNPs were filtered out
SNPS_FILTERED=$((TOTAL_SNPS - FILTERED_SNPS))
echo "Number of SNPs filtered out: $SNPS_FILTERED"
