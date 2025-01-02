#!/bin/bash 

input_dir="mounted-data-readonly/ms_gwas_results"    
out_dir="gwas_results"

pheno_list=("primary_progressive_ms_coded" "ms_all_comers" "primary_progressive_vs_rrsms_and_spms" "relapsing_remitting_ms_coded" "rrms_with_spms" "secondary_progressive_ms_coded")

# Loop over each trait
for trait in "${pheno_list[@]}"; do
    echo "Combining results for $trait"

    # Create an output file for the combined results
    output_file="${out_dir}/SPACox_${trait}.cox"

    # Remove the output file if it already exists
    rm -f "$output_file"

    # Combine the Regenie result files for all 22 chromosomes
    for chr in {1..22}; do
        input_file="${input_dir}/SPACox_${trait}_chr${chr}_results.txt"
        
        # Append the content of the chromosome file to the output file
        # Assuming first file contains header, we skip header from other files
        if [ "$chr" -eq 1 ]; then
            cat "$input_file" >> "$output_file"
        else
            tail -n +2 "$input_file" >> "$output_file"
        fi
    done

    echo "Results combined for $trait into $output_file"
done

