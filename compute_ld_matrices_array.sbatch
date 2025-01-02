#!/bin/bash
#SBATCH --job-name=compute_ld_matrices
#SBATCH --output=logs/compute_ld_matrix_%A_%a.out
#SBATCH --error=logs/compute_ld_matrix_%A_%a.err
#SBATCH --time=4:00:00              # Adjust time limit as needed
#SBATCH --cpus-per-task=4           # Number of CPUs per task
#SBATCH --mem=16G                   # Memory per task
#SBATCH --partition=compute         # Adjust partition/queue name as per your cluster

# Load necessary modules
#module load plink/2.0
#module load plink/1.9
#module load bedtools

# Get the block ID for this task
BLOCK_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ld_blocks/block_ids.txt)

# Run the compute LD matrix script for this block
./compute_ld_matrix_single_block.sh "$BLOCK_ID"

