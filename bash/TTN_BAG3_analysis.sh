#!/bin/bash

# To calculate maf, use sjlife.fam on /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3. Use pheno/sjlife_ttn_bag3.pheno
awk '{print $1"\t"$2}' pheno/sjlife_ttn_bag3.pheno > samples_for_maf.txt

module load plink/1.90b
plink --bfile sjlife --keep samples_for_maf.txt --freq --out sjlife.freq.out

# First removing trailing spaces in the file, then removing chr
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' sjlife.freq.out.frq | awk -F " " '{gsub(/chr/, "", $2); print}' OFS="\t" > sjlife.freq.out.frq_edited1

awk -F'\t' '{print $0 FS $2":"$4":"$3; exit}' Summary_results_AN.txt
awk -F'\t' '{print $0"\t"$2}' Summary_results_AN.txt