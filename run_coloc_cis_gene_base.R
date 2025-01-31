#-------------------------------------------------- 
# 0. Load packages
#--------------------------------------------------
library(coloc)
library(data.table)

#--------------------------------------------------
# 1. Read in your (gzipped) summary stats as data.tables
#--------------------------------------------------
gwas_data <- fread("/mnt/project/genuity_data/time2event_results/ms_rrms_age_onset.regenie.gz")
eqtl_data <- fread("/mnt/project/publically_available_supporting_files/rosmap_brain/celltype-eqtl-sumstats.Ast.tsv.gz")

#--------------------------------------------------
# 2. Format eQTL SNP ID to chr:pos:REF:ALT
#--------------------------------------------------
eqtl_data[, chr38 := gsub("chr", "", chr38)]  # Remove "chr" prefix
eqtl_data[, chr38 := as.numeric(chr38)]
etal_data <- eqtl_data[!is.na(chr38)] # remove rows for X,Y, and MT
eqtl_data[, formatted_snp := paste(chr38, pos38, REF, ALT, sep = ":")]

#--------------------------------------------------
# 3. List of unique genes
#--------------------------------------------------
all_genes <- unique(eqtl_data$gene_symbol)

#--------------------------------------------------
# 4. Initialize a list to store results
#--------------------------------------------------
coloc_results_list <- list()  

# Sample sizes (replace with actual values)
GWAS_N  <- 26119  # GWAS total sample size
eQTL_N  <-  576   # eQTL total sample size

#--------------------------------------------------
# 5. Loop over each gene
#--------------------------------------------------
for (gene in all_genes) {
  
  # Subset eQTL data for this gene
  eqtl_sub <- eqtl_data[gene_symbol == gene]
  
  # If no eQTL SNPs, skip
  if (nrow(eqtl_sub) < 1) next
  
  # Identify chromosome
  chr_of_interest <- unique(eqtl_sub$chr38)
  if (length(chr_of_interest) != 1) next  # Skip if multiple chromosomes for a gene
  
  # Subset GWAS to matching chromosome and positions
  eqtl_positions <- unique(eqtl_sub$pos38)
  gwas_sub <- gwas_data[CHROM == chr_of_interest & GENPOS %in% eqtl_positions]

  # If no overlapping GWAS SNPs, skip
  if (nrow(gwas_sub) < 1) next
  
  #--------------------------------------------------
  # 6. Align SNPs & Check Strand Flipping
  #--------------------------------------------------
  
  # Merge GWAS and eQTL by chromosome & position
  merged_snps <- merge(eqtl_sub, gwas_sub, 
                       by.x = c("chr38", "pos38"), 
                       by.y = c("CHROM", "GENPOS"))

  if (nrow(merged_snps) < 1) next  # Skip if no shared SNPs
  
  # Detect strand flips: when REF/ALT are swapped
  merged_snps[, strand_flipped := (REF == ALLELE1 & ALT == ALLELE0)]

  # Apply strand flipping correction by negating beta
  merged_snps[strand_flipped == TRUE, `:=` (
      beta= -1 * beta,
      ALT_AF = 1-ALT_AF,
      formatted_snp = paste(chr38, pos38, ALT, REF, sep = ":")
  )]
 merged_snps <- unique(merged_snps, by = "formatted_snp")
    
  #--------------------------------------------------
  # 7. Prepare Data for Coloc
  #--------------------------------------------------
  eqtl_dataset <- list(
    snp     = merged_snps$formatted_snp, 
    beta    = merged_snps$beta,  
    varbeta = (merged_snps$se)^2, 
    MAF     = merged_snps$ALT_AF,  
    N       = eQTL_N,  
    type    = "quant"
  )
  
  gwas_dataset <- list(
    snp     = merged_snps$ID,  
    beta    = merged_snps$BETA,  
    varbeta = (merged_snps$SE)^2,  
    MAF     = merged_snps$A1FREQ,  
    N       = GWAS_N,  
    type    = "quant"
  )
  
  #--------------------------------------------------
  # 8. Run Coloc
  #--------------------------------------------------
 suppressWarnings(suppressMessages({
  this_coloc <- coloc.abf(eqtl_dataset, gwas_dataset)
  }))
  
  # Extract results
  summary_vec <- this_coloc$summary
  
  # Store results
  coloc_results_list[[gene]] <- data.table(
    gene_symbol = gene,
    PP.H0 = summary_vec["PP.H0.abf"],
    PP.H1 = summary_vec["PP.H1.abf"],
    PP.H2 = summary_vec["PP.H2.abf"],
    PP.H3 = summary_vec["PP.H3.abf"],
    PP.H4 = summary_vec["PP.H4.abf"]
  )
}

#--------------------------------------------------
# 9. Combine all gene results into a single data.table
#--------------------------------------------------
coloc_results <- rbindlist(coloc_results_list)

# Sort by highest colocalization probability (PP.H4)
#setorder(coloc_results, -PP.H4.abf)

#--------------------------------------------------
# 10. Save or export results
#--------------------------------------------------
fwrite(coloc_results, "coloc_ms_rrms_age_onset_vs_Ast.txt", sep = "\t", quote = FALSE)

# View top colocalized genes
head(coloc_results)
