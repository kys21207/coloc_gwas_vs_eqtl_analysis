# run_coloc_susie_parallel.R 
# R script to run coloc.susie for a given block using parallel execution

# Load required libraries
suppressMessages(library(coloc))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(stringr))
suppressMessages(library(tidyverse))
suppressMessages(library(doMC))
suppressMessages(library(tictoc))
suppressMessages(library(Rfast))

# Input arguments
args <- commandArgs(trailingOnly = TRUE)

block_id   <- args[1]
gwasA_file <- args[2]
gwasB_path <- args[3]
ld_matrix_file <- args[4]
output_file <- args[5]

tic("Coloc execution")
# Load gwasB files from coloc_input folder
pattern <- paste0("ALL-GWAS_ALL-GWAS_.*\\.regenie_block_", block_id, "\\.txt")
files <- list.files(path = gwasB_path, pattern = pattern , full.names = TRUE)
#file_names <- basename(files)
#cat(gwasB_path, "\n")
#cat(block_id,"\n")
cat("number of files", length(files), "\n")

# Load input data
LD_MATRICES_DIR="filesystems/ld_matrices"
dataset1 <- read.table(gwasA_file, header=TRUE)
ld_matrix <- as.matrix(fread(ld_matrix_file, header=FALSE))

# Load SNP IDs from BIM file
bim_file <- file.path(LD_MATRICES_DIR, paste0("ld_matrix_block_", block_id, ".bim"))
bim <- fread(bim_file, header = FALSE)
snp_ids <- bim$V2

# Register parallel backend
#n_cores <- detectCores() - 1 # Use one less core than available
n_cores <- 16
cat("number of cores", n_cores, "\n")
registerDoMC(n_cores)

# Run coloc analysis in parallel
results <- foreach(gwasB_file = files, .combine = rbind, .packages = c("data.table", "coloc", "stringr", "tidyverse")) %dopar% {
  # Load dataset2
  dataset2 <- read.table(gwasB_file, header=TRUE)
 # cat("file names:",basename(gwasB_file),"\n")
    
  # Identify common SNPIDs across all datasets
  common_SNPIDs <- Reduce(intersect, list(dataset1$ID, dataset2$ID, snp_ids))

  # Filter datasets to include only common SNPIDs
  dataset1_filtered <- dataset1[match(common_SNPIDs, dataset1$ID), ]
  dataset2_filtered <- dataset2[match(common_SNPIDs, dataset2$ID), ]

  ld_matrix_filtered <- ld_matrix[match(common_SNPIDs, snp_ids), match(common_SNPIDs, snp_ids)]
  colnames(ld_matrix_filtered) <- common_SNPIDs
  rownames(ld_matrix_filtered) <- common_SNPIDs

  # Prepare datasets for coloc.susie
  coloc_dataset1 <- list(beta=dataset1_filtered$BETA, varbeta=dataset1_filtered$SE^2, position=dataset1_filtered$GENPOS, snp=dataset1_filtered$ID, type="cc", N=dataset1_filtered$N[1], LD=ld_matrix_filtered)
  coloc_dataset2 <- list(beta=dataset2_filtered$BETA, varbeta=dataset2_filtered$SE^2, position=dataset2_filtered$GENPOS, snp=dataset2_filtered$ID, type="cc", N=dataset2_filtered$N[1], LD=ld_matrix_filtered)
suppressMessages({
  suppressWarnings({
  # Run coloc.abf and coloc.susie
  result <- coloc.abf(coloc_dataset1, coloc_dataset2)
  extracted_part <- str_extract(gwasB_file, "ENSG[0-9]+")
  result <- as.data.frame(t(result$summary)) %>% mutate(gene=extracted_part, causal="single", block=block_id, .before=nsnps)
  if(result$PP.H4.abf > 0.7) {
      result1 <- coloc.susie(coloc_dataset1, coloc_dataset2)
      if (!is.null(result1$summary)) {
            result <- as.data.frame(t(result1$summary)) %>% mutate(gene=extracted_part, causal="multi", block=block_id, .before=nsnps)
      }
    }
      })
    })

  return(result)
}

# Save results
if (file.exists(output_file)) {
  write.table(results, paste(output_file), quote=F, row.names=F, sep="\t", append=TRUE, col.names=F)
} else {
  write.table(results, paste(output_file), quote=F, row.names=F, sep="\t")
}

cat("Completed coloc.susie for", gwasA_file, "and eQTL\n")
toc()

