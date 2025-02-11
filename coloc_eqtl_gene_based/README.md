## 1. bedtools 
It will be very useful to intersect between two bed files. <br>
For example, pull out a subset from a gwas/eqtl sumstats according to LD block regions. <br>
## 2. coloc analysis workflow
1. select LD blocks containing the significant GWAS snps that are filtered by the pre-specified threshold (filter_gwas+sumstats.sh & identify_LD_blocks_based_sig_gwas.sh) 
2. Pull a subset eqtl sumstats out according to the selected LD blocks (filter_eqtl_sumstats.sh)
3. run coloc based on a gwas and a subset of eqtl sumstats (coloc_gene_base_onek1k.r)  
