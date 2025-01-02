#!/bin/bash 

#chro=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
chro=11 

for chr in "${chro[@]}"; do
    plink2 --vcf kgp/kgp_renamed_chrom${chr}.vcf.gz --remove kgp/non_eur_ids.txt --make-pgen --out kgp/kgp_chr${chr}
done

