#!/bin/bash

# Check if arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Wrong number of arguments supplied."
    echo "Usage: $0 gwas_dir eqtl_dir"
    echo "Example: $0 /mnt/project/harmonized_genuity_results /mnt/project/publically_available_supporting_files/rosmap_brain"
    exit 1
fi

# Input directories from command line arguments
GWAS_INPUT_DIR=$1
EQTL_INPUT_DIR=$2
TIMESTAMP="2025-02-07_15-16-14"  # Updated to current time
CURRENT_USER="kys21207"
PROJECT="project-Gv45qjQ09Vk2p6X7q5xJ42PV"

# Log input parameters
echo "=== Analysis Parameters ==="
echo "GWAS Directory: $GWAS_INPUT_DIR"
echo "eQTL Directory: $EQTL_INPUT_DIR"
echo "Current User: $CURRENT_USER"
echo "Timestamp: $TIMESTAMP"
echo "Project: $PROJECT"
echo "==========================="

# Process one eQTL file at a time
for eqtl_file in $EQTL_INPUT_DIR/*.all.tsv.gz; do
    eqtl_name_prefix=$(basename "$eqtl_file" .all.tsv.gz)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing eQTL file: $eqtl_name_prefix"

    # Create batch directory for this eQTL
    BATCH_DIR="${PROJECT}:/analysis_KJ/coloc_gwas_eqtls/results"
    
    # Remove "/mnt/project/" from the paths
    eqtl_cleaned_path=${eqtl_file#/mnt/project/}

    # Create command for processing all GWAS files
    CMD="set -e; mkdir -p /data /results && "
    
    # Add copy commands for input files
    CMD+="cp /mnt/project/publically_available_supporting_files/onek1k/oneK1K_samplesize.csv /data/ && "
    CMD+="cp /mnt/project/${eqtl_cleaned_path} /data/ && "
    
    # Add commands for each GWAS file
    for gwas_file in $GWAS_INPUT_DIR/*.regenie.tsv.gz; do
        gwas_name_prefix=$(basename "$gwas_file" .regenie.tsv.gz)
        gwas_cleaned_path=${gwas_file#/mnt/project/}
        CMD+="cp /mnt/project/${gwas_cleaned_path} /data/ && "
        CMD+="cd /data && Rscript /opt/notebooks/codes/run_coloc_gene_base.r ${gwas_name_prefix}.regenie.tsv.gz ${eqtl_name_prefix}.all.tsv.gz oneK1K_samplesize.csv && "
    done

    # Remove the trailing " && "
    CMD=${CMD% && }

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting batch job for all GWAS files vs $eqtl_name_prefix"

    # Submit swiss-army-knife job
    dx run swiss-army-knife \
        -iin="${PROJECT}:/publically_available_supporting_files/onek1k/oneK1K_samplesize.csv" \
        -iin="${PROJECT}:/${eqtl_cleaned_path}" \
        -iimage_file="${PROJECT}:/analysis_KJ/coloc_gwas_eqtls/code/coloc-analysis.tar" \
        -icmd="${CMD}" \
        --instance-type mem3_ssd1_v2_x8 \
        --destination="${BATCH_DIR}" \
        --brief \
        --yes

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed submission for eQTL: $eqtl_name_prefix"
    echo "----------------------------------------"
done

echo "All coloc analyses submitted."
