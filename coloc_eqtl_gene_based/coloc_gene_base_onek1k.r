#-------------------------------------------------- 
# 0. Load packages
#--------------------------------------------------
library(coloc)
library(data.table)

#--------------------------------------------------
# 1. Read in your (gzipped) summary stats as data.tables
#--------------------------------------------------

# Input arguments
args <- commandArgs(trailingOnly = TRUE)

eqtl_file        <- args[1]
gwas_dir         <- args[2]
output_dir       <- args[3]

eqtl_recode <- read.table("/mnt/project/publically_available_supporting_files/onek1k/oneK1K_samplesize.csv",header=T,sep=",")
#gwas_files <-"/mnt/project/genuity_data/time2event_results/ms_all_comers_age_onset.regenie.gz"
gwas_files <- list.files(path = gwas_dir, pattern = "\\.regenie.tsv.gz$", full.names = TRUE)
# oneK1K data parameters
#molecular_trait_id      chromosome      position        ref     alt     variant ma_samples      maf     pvalue  beta    se      type    ac      an      r2      molecular_trait_object_id     gene_id median_tpm      rsid
eqtl_file <- "/mnt/project/publically_available_supporting_files/onek1k/cDC2.all.tsv.gz"
eqtl_data <- fread(eqtl_file)
eqtl_file_prefix <- sub("\\.all\\.tsv\\.gz$", "", basename(eqtl_file))
#eqtl_cell <- sub("^celltype-eqtl-sumstats\\.", "", eqtl_file_prefix)
eqtl_cell <-eqtl_file_prefix

eqtl_sample_size <- subset(eqtl_recode, sample_group == eqtl_file_prefix)

#--------------------------------------------------
# 2. Format eQTL SNP ID to chr:pos:REF:ALT
#--------------------------------------------------
#eqtl_data[, chr38 := gsub("chr", "", chr38)]  # Remove "chr" prefix
eqtl_data[, chromosome := as.numeric(chromosome)]
eqtl_data <- eqtl_data[!is.na(chromosome)] # remove rows for X,Y, and MT
eqtl_data[, formatted_snp := paste(chromosome, position, ref, alt, sep = ":")]

#--------------------------------------------------
# 3. List of unique genes
#--------------------------------------------------
all_genes <- unique(eqtl_data$molecular_trait_id)

for (gwas in gwas_files){
print(gwas)
gwas_data <- fread(gwas)

#--------------------------------------------------
# 4. Initialize a list to store results
#--------------------------------------------------
coloc_results_list <- list()  

# Sample sizes (replace with actual values)
GWAS_N  <- gwas_data$N[1]  # GWAS total sample size
eQTL_N  <-  eqtl_sample_size$sample_size   # eQTL total sample size

#--------------------------------------------------
# 5. Loop over each gene
#--------------------------------------------------
for (gene in all_genes) {
  
  # Subset eQTL data for this gene
  eqtl_sub <- eqtl_data[molecular_trait_id == gene]
  
  # If no eQTL SNPs, skip
  if (nrow(eqtl_sub) < 1) next
  
  # Identify chromosome
  chr_of_interest <- unique(eqtl_sub$chromosome)
  if (length(chr_of_interest) != 1) next  # Skip if multiple chromosomes for a gene
  
  # Subset GWAS to matching chromosome and positions
  eqtl_positions <- unique(eqtl_sub$position)
  gwas_sub <- gwas_data[CHROM == chr_of_interest & GENPOS %in% eqtl_positions]

  # If no overlapping GWAS SNPs, skip
  if (nrow(gwas_sub) < 1) next
  
  #--------------------------------------------------
  # 6. Align SNPs & Check Strand Flipping
  #--------------------------------------------------
  
  # Merge GWAS and eQTL by chromosome & position
  merged_snps <- merge(eqtl_sub, gwas_sub, 
                       by.x = c("chromosome", "position"), 
                       by.y = c("CHROM", "GENPOS"))

  if (nrow(merged_snps) < 1) next  # Skip if no shared SNPs
  
  # Detect strand flips: when REF/ALT are swapped
  merged_snps[, strand_flipped := (ref == ALLELE1 & alt == ALLELE0)]

  # Apply strand flipping correction by negating beta
  merged_snps[strand_flipped == TRUE, `:=` (
      beta= -1 * beta,
      maf = 1-maf,
      formatted_snp = paste(chromosome, position, alt, ref, sep = ":")
  )]
 merged_snps <- unique(merged_snps, by = "formatted_snp")
 merged_snps <- unique(merged_snps, by = "ID")
    
  #--------------------------------------------------
  # 7. Prepare Data for Coloc
  #--------------------------------------------------
  eqtl_dataset <- list(
    snp     = merged_snps$formatted_snp, 
    beta    = merged_snps$beta,  
    varbeta = (merged_snps$se)^2, 
    MAF     = merged_snps$maf,  
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
  capture.output({
    this_coloc <- coloc.abf(eqtl_dataset, gwas_dataset)
  }, file = "/dev/null")
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
gwas_file_prefix <- sub("\\.regenie(\\.tsv)?\\.gz$", "", basename(gwas))    

fwrite(coloc_results, paste0(output_dir,"/coloc_",gwas_file_prefix,"_vs_",eqtl_cell,".txt"), sep = "\t", quote = FALSE)

# View top colocalized genes
#head(coloc_results)
}
