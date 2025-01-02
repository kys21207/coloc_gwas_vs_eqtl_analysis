# run_significant_blocks.sh
# Script to run coloc.susie for each significant LD block

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a file name for the final results."
    echo "Usage: $0 file_name"
    exit 1
fi

finaloutputname=$1
batch=$2

# Input Arguments
COLOC_INPUT_DIR="coloc_inputs/${batch}" # Directory containing coloc input files
RESULTS_DIR="coloc_results/${batch}" # Directory where coloc results will be stored
gwasB_path=$COLOC_INPUT_DIR

mkdir -p $RESULTS_DIR

# Load R environment if needed (e.g., conda activate or module load R)
# source activate R_env

# Loop through each significant block and run coloc.susie for GWAS A vs all GWAS B
for block_id_file in $COLOC_INPUT_DIR/gwasA_block_*.txt; do
  block_id=$(basename $block_id_file | sed 's/gwasA_block_\(.*\)\.txt/\1/')
  
    Rscript codes/run_coloc_susie.R "$block_id" "$block_id_file" "$gwasB_path" "filesystems/ld_matrices/ld_matrix_block_${block_id}.ld.gz" "$RESULTS_DIR/coloc_result_${finaloutputname}_vs_ALL_eQTL_${batch}.txt"

    echo "Completed coloc.susie for GWAS A, GWAS B Block ID: $block_id"
done

echo "All coloc.susie analyses completed."
rm $COLOC_INPUT_DIR/gwasA_block_*.txt

