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



# (base) [aneupane@splprhpc09 ttn_bag3]$ head sjlife_results.assoc.logistic.clean.Psorted
# CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P
# 10 chr10:119570337:A:C 119570337 C ADD 1645 4.403 0.4072 1.982 9.782 3.64 0.0002725
# 10 chr10:119640101:T:C 119640101 C ADD 1630 3.383 0.348 1.71 6.691 3.502 0.0004613
# 10 chr10:119644280:G:A 119644280 A ADD 1645 3.318 0.3462 1.683 6.539 3.463 0.0005332
# 10 chr10:119754685:C:G 119754685 G ADD 1645 3.166 0.3497 1.595 6.283 3.295 0.0009831
# 10 chr10:119833210:A:G 119833210 G ADD 1645 2.973 0.3447 1.513 5.843 3.161 0.001571


# (base) [aneupane@splprhpc09 ttn_bag3]$ head sjlife_results.assoc.logistic.clean.Psorted.formetal
# CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA
# 10 chr10:119570337:A:C 119570337 C ADD 1645 4.403 0.4072 1.982 9.782 3.64 0.0002725 A C A 1.48229
# 10 chr10:119640101:T:C 119640101 C ADD 1630 3.383 0.348 1.71 6.691 3.502 0.0004613 T C T 1.21876
# 10 chr10:119644280:G:A 119644280 A ADD 1645 3.318 0.3462 1.683 6.539 3.463 0.0005332 G A G 1.19936


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


# (base) [aneupane@splprhpc09 GWAS]$ awk 'NR==FNR{a[$1":"$3]=$0;next}($1":"$3 in a){print $0, a[$1":"$3]}' DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal \
# > DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal| head
# CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA
# 19 chr19:44908822:C:T 44908822 T ADD 494 0.570581  0.09213 -0.7417 -0.3805 -6.09 2.326e-09 C T C -0.5611 19 chr19:44908822:C:T 44908822 T ADD 2949 0.713124  0.04545 -0.4272 -0.2491 -7.44 1.314e-13 C T C -0.3381
# 19 chr19:44897490:T:A 44897490 A ADD 495 0.595353  0.09386 -0.7026 -0.3346 -5.525 5.431e-08 T A T -0.518601 19 chr19:44897490:T:A 44897490 A ADD 2956 0.796045  0.07043 -0.3661 -0.09002 -3.238 0.001217 T A T -0.2281
# 19 chr19:44899005:T:G 44899005 G ADD 494 0.60236  0.09211 -0.6874 -0.3263 -5.503 6.114e-08 T G T -0.5069 19 chr19:44899005:T:G 44899005 G ADD 2952 0.819878  0.06985 -0.3355 -0.06168 -2.843 0.0045 T G T -0.1986
# 19 chr19:44909976:G:T 44909976 T ADD 490 0.614098  0.09077 -0.6655 -0.3096 -5.371 1.231e-07 G T G -0.487601 19 chr19:44909976:G:T 44909976 T ADD 2953 0.709354  0.04503 -0.4317 -0.2552 -7.626 3.249e-14 G T G -0.343401


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


# (base) [aneupane@splprhpc09 ttn_bag3]$ head ccss_exp_results.assoc.logistic.clean.Psorted.formetal.common
# CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA
# 10 10:119684833 119684833 T ADD 1483 5.115 0.3806 2.426 10.79 4.289 1.798e-05 G T G 1.63218
# 10 10:119800012 119800012 T ADD 1449 3.241 0.2758 1.888 5.564 4.263 2.013e-05 C T C 1.17588
# 10 10:119794262 119794262 T ADD 1483 3.192 0.2729 1.87 5.449 4.253 2.106e-05 C T C 1.16065
# 10 10:119875098 119875098 C ADD 1484 3.457 0.337 1.786 6.691 3.681 0.0002328 G C G 1.2404
# 10 10:119907508 119907508 C ADD 1484 1.891 0.2027 1.271 2.813 3.142 0.001678 T C T 0.637106
# 10 10:119886263 119886263 A ADD 1483 1.867 0.2124 1.231 2.831 2.94 0.003278 G A G 0.624333
# 10 10:119877120 119877120 T ADD 1482 1.866 0.2124 1.231 2.83 2.938 0.003301 C T C 0.623797
# 10 10:119882821 119882821 T ADD 1482 1.866 0.2124 1.231 2.829 2.937 0.003318 C T C 0.623797
# 10 10:119877662 119877662 CA ADD 1484 1.862 0.2121 1.229 2.822 2.932 0.003365 C CA C 0.621651



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
