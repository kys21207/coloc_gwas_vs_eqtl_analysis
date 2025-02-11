#!/bin/bash

# Ensure the input directory and output directory are provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 input_directory output_directory"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR="/opt/notebooks/data"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Use tmp as the temporary directory for intermediate files
TMP_DIR=tmp

# Create the temporary directory if it doesn't exist
mkdir -p $TMP_DIR

# Function to convert a single .tsv.gz file to .bed.gz
convert_to_bed() {
    local input_file=$1
    local output_file=$2
    zcat $input_file | tail -n +2 | awk 'BEGIN {OFS="\t"} {print $2, $3, $3, $1, $2, $3, $4, $5, $8, $10, $11}' | gzip > $output_file
}

# Find all .all.tsv.gz files and convert them to .bed.gz
find $INPUT_DIR -name "*.all.tsv.gz" | while read -r file; do
    base_name=$(basename "$file" .all.tsv.gz)
    output_file="$OUTPUT_DIR/${base_name}.bed.gz"
    convert_to_bed "$file" "$output_file"
done

# Clean up temporary files
rm -r $TMP_DIR

echo "All .all.tsv.gz files have been converted to .bed.gz"
