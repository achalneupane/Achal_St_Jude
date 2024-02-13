#!/bin/bash

EUR=DYSLPDM_chltotl_rank_chr*.assoc.linear


## Process PLINK result files for meta-analysis using METAL; add A2 and BETA
# SJLIFE
# Find the header file and save it to a variable
header_file=$(ls DYSLPDM_chltotl_rank_chr*.assoc.linear | head -n 1)
# Concatenate all files, skipping the header for each file except the first
cat $header_file > all_chr_EUR  # Copy the header from the first file
cat $(ls DYSLPDM_chltotl_rank_chr*.assoc.linear | sort -V) | grep -v CHR >> all_chr_EUR

tr '\t' ' ' < all_chr_EUR | sed 's/  */ /g' > DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean
rm all_chr_EUR
{ head -n 1 DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean && tail -n +2 DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean | sort -k12,12g; } > DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted


AFR=DYSLPDM_chltotl_rankeur_afr_chr*.assoc.linear

header_file=$(ls DYSLPDM_chltotl_rankeur_afr_chr*.assoc.linear | head -n 1)
# Concatenate all files, skipping the header for each file except the first
cat $header_file > all_chr_AFR  # Copy the header from the first file
cat $(ls DYSLPDM_chltotl_rankeur_afr_chr*.assoc.linear | sort -V) | grep -v CHR >> all_chr_AFR

tr '\t' ' ' < all_chr_AFR | sed 's/  */ /g' > DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean
rm all_chr_AFR
{ head -n 1 DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean && tail -n +2 DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean | sort -k12,12g; } > DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted




awk 'BEGIN {OFS="\t"} NR==1 {print $0, "OR"} NR>1 {print $0, exp($7)}' DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted > tmp
awk 'BEGIN {OFS=" "} {print $1, $2, $3, $4, $5, $6, $13, $14, $8, $9, $10, $11, $12}' tmp > DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted

awk 'BEGIN {OFS="\t"} NR==1 {print $0, "OR"} NR>1 {print $0, exp($7)}' DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted > tmp
awk 'BEGIN {OFS=" "} {print $1, $2, $3, $4, $5, $6, $13, $14, $8, $9, $10, $11, $12}' tmp > DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted
rm tmp

# EUR
awk '{split($2, a, ":"); print $0, a[3], a[4]}' DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted \
| awk '{if($4==$13) A2=$14; else A2=$13; BETA=log($7); print $0, A2, BETA}' \
| awk 'BEGIN{print "CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA"}FNR>1{print}' \
> DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal
# AFR
awk '{split($2, a, ":"); print $0, a[3], a[4]}' DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted \
| awk '{if($4==$13) A2=$14; else A2=$13; BETA=log($7); print $0, A2, BETA}' \
| awk 'BEGIN{print "CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA"}FNR>1{print}' \
> DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal



## Identify variants that match between the 2 datasets, based on chr:pos and then based on both A1 and A2
awk 'NR==FNR{a[$1":"$3]=$0;next}($1":"$3 in a){print $0, a[$1":"$3]}' DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal \
DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal| awk '($4==$20 || $4==$31) && ($15==$20 || $15==$31)' > common_variants_among_2_datasets.txt

## Fix space
sed 's/[[:blank:]]\+/ /g' common_variants_among_2_datasets.txt > common_variants_among_2_datasets_updated.txt

# Then prepare final files for each dataset for METAL
awk '!a[$1":"$3]++' common_variants_among_2_datasets_updated.txt | cut -d' ' -f1-16 | awk '{if(FNR>1) $2=$1":"$3; print}' \
> DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal.common
awk '!a[$1":"$3]++' common_variants_among_2_datasets_updated.txt | cut -d' ' -f17-32 | awk '{if(FNR>1) $2=$1":"$3; print}' \
> DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal.common



## Prepare the config file for METAL
# run_metal.script

## Run METAL
~/bin/generic-metal/metal run_metal.script
SCHEME STDERR
STDERR SE

MARKER SNP
ALLELE A1 A2
EFFECT BETA
PVALUE P
PROCESS DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal.common

PROCESS DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal.common

OUTFILE Meta_DYSLPDM_chltotl_EUR_AFR_fixed_ .tbl

ANALYZE HETEROGENEITY

QUIT



## A total of 7319447 variants present in both datasets were meta-analyzed; process the results further
