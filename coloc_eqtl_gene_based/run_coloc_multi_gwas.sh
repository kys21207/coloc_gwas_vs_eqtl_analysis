#!/bin/bash

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a gwas directory and a eqtl directory."
    echo "Usage: $0 file_name"
    exit 1
fi

# examples
#EQTL_INPUT_DIR="/mnt/project/publically_available_supporting_files/rosmap_brain"
#GWAS_INPUT_DIR="/mnt/project/harmonized_genuity_results"

GWAS_INPUT_DIR=$1
EQTL_INPUT_DIR=$2
EQTL_INPUT_SAMPLE_SIZE_FILE=$3

# Input Arguments
RESULTS_DIR="/opt/notebooks/coloc_results" # Directory where coloc results will be stored

mkdir -p $RESULTS_DIR

chmod +x ./codes/filter_gwas_sumstats.sh
chmod +x ./codes/identify_sig_LD_blockID.sh
chmod +x ./codes/filter_eqtl_sumstats.sh

# Load R environment if needed (e.g., conda activate or module load R)
# source activate R_env

for gwas_file in $GWAS_INPUT_DIR/*.regenie.tsv.gz; do
  gwas_name_prefix=$(basename "$gwas_file" .regenie.tsv.gz)   
  
  ./codes/filter_gwas_sumstats.sh $gwas_file $gwas_name_prefix
  ./codes/identify_sig_LD_blockID.sh $gwas_name_prefix
# Loop through each significant block and run coloc.susie for GWAS A vs all GWAS B
for eqtl_file in $EQTL_INPUT_DIR/*.bed.gz; do
   eqtl_name_prefix=$(basename "$eqtl_file" .bed.gz) 

   ./codes/filter_eqtl_sumstats.sh /opt/notebooks/tmp_output/significant_ld_blocks.txt $eqtl_file $eqtl_name_prefix
  if [ -f "$eqtl_file" ]; then
    # Process the file with fread in R
    Rscript ./codes/run_coloc_gene_base.r "${gwas_file}" "/opt/notebooks/tmp_input/${eqtl_name_prefix}_subset.tsv.gz" "$EQTL_INPUT_SAMPLE_SIZE_FILE" "$RESULTS_DIR"
    
#    dx upload /opt/notebooks/coloc_results/*.txt --destination project-Gv45qjQ09Vk2p6X7q5xJ42PV:/analysis_KJ/coloc_gwas_eqtls/results/
    rm /opt/notebooks/coloc_results/*.txt

  else
    echo "File '$eqtl_file' does not exist or is non-readable."
  fi 
 #   echo "Completed coloc for GWAS $gwas_name_prefix, eQTLs"
    echo " gwas name: $eqtl_file "
done
rm /opt/notebooks/tmp_input/*.gz
done

echo "All coloc analyses completed."
#rm $COLOC_INPUT_DIR/gwasA_block_*.tx

echo "All coloc analyses completed."
#rm $COLOC_INPUT_DIR/gwasA_block_*.txt

