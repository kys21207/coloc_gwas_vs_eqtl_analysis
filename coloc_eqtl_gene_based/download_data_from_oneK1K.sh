#!/bin/bash

base_url="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038"
for i in $(seq -w 606 629); do
    wget "${base_url}/QTD000${i}/QTD000${i}.all.tsv.gz"
done

declare -A rename_map=(
    [QTD000606]="B_intermediate"
    [QTD000607]="B_memory"
    [QTD000608]="B_naive"
    [QTD000609]="CD14_Mono"
    [QTD000610]="CD16_Mono"
    [QTD000611]="CD4_CTL"
    [QTD000612]="CD4_Naive"
    [QTD000613]="CD4_TCM"
    [QTD000614]="CD4_TEM"
    [QTD000615]="CD8_Naive"
    [QTD000616]="CD8_TCM"
    [QTD000617]="CD8_TEM"
    [QTD000618]="HSPC"
    [QTD000619]="MAIT"
    [QTD000620]="NK"
    [QTD000621]="NK_CD56bright"
    [QTD000622]="NK_Proliferating"
    [QTD000623]="Plasmablast"
    [QTD000624]="Platelet"
    [QTD000625]="Treg"
    [QTD000626]="cDC2"
    [QTD000627]="dnT"
    [QTD000628]="gdT"
    [QTD000629]="pDC"
)

for dataset_id in "${!rename_map[@]}"; do
    mv "${dataset_id}.all.tsv.gz" "${rename_map[$dataset_id]}.all.tsv.gz"
done
