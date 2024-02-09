#!/bin/bash

## Process PLINK result files for meta-analysis using METAL; add A2 and BETA
# (base) [aneupane@splprhpc09 ttn_bag3]$ head sjlife_results.assoc.logistic.clean.Psorted
# CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P
# 10 chr10:119570337:A:C 119570337 C ADD 1645 4.403 0.4072 1.982 9.782 3.64 0.0002725
# 10 chr10:119640101:T:C 119640101 C ADD 1630 3.383 0.348 1.71 6.691 3.502 0.0004613
# 10 chr10:119644280:G:A 119644280 A ADD 1645 3.318 0.3462 1.683 6.539 3.463 0.0005332
# 10 chr10:119754685:C:G 119754685 G ADD 1645 3.166 0.3497 1.595 6.283 3.295 0.0009831
# 10 chr10:119833210:A:G 119833210 G ADD 1645 2.973 0.3447 1.513 5.843 3.161 0.001571
# 10 chr10:119691120:G:C 119691120 C ADD 1644 1.426 0.1185 1.13 1.798 2.992 0.002776
# 2 chr2:178395655:A:G 178395655 A ADD 1643 0.7096 0.1157 0.5657 0.8902 -2.966 0.003017
# 10 chr10:119610713:G:A 119610713 A ADD 1645 2.82 0.3518 1.415 5.619 2.947 0.003205
# 10 chr10:119877390:G:A 119877390 A ADD 1645 2.966 0.3709 1.434 6.136 2.932 0.003371

tr '\t' ' ' < sjlife_results.assoc.logistic | sed 's/  */ /g' > sjlife_results.assoc.logistic.clean
{ head -n 1 sjlife_results.assoc.logistic.clean.Psorted && tail -n +2 sjlife_results.assoc.logistic.clean.Psorted | sort -k12,12n; } > sjlife_results.assoc.logistic.clean.Psorted.sorted


# SJLIFE
awk '{split($2, a, ":"); print $0, a[3], a[4]}' sjlife_results.assoc.logistic.clean.Psorted \
| awk '{if($4==$13) A2=$14; else A2=$13; BETA=log($7); print $0, A2, BETA}' \
| awk 'BEGIN{print "CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA"}FNR>1{print}' \
> sjlife_results.assoc.logistic.clean.Psorted.formetal

# (base) [aneupane@splprhpc09 ttn_bag3]$ head sjlife_results.assoc.logistic.clean.Psorted.formetal
# CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA
# 10 chr10:119570337:A:C 119570337 C ADD 1645 4.403 0.4072 1.982 9.782 3.64 0.0002725 A C A 1.48229
# 10 chr10:119640101:T:C 119640101 C ADD 1630 3.383 0.348 1.71 6.691 3.502 0.0004613 T C T 1.21876
# 10 chr10:119644280:G:A 119644280 A ADD 1645 3.318 0.3462 1.683 6.539 3.463 0.0005332 G A G 1.19936
# 10 chr10:119754685:C:G 119754685 G ADD 1645 3.166 0.3497 1.595 6.283 3.295 0.0009831 C G C 1.15247
# 10 chr10:119833210:A:G 119833210 G ADD 1645 2.973 0.3447 1.513 5.843 3.161 0.001571 A G A 1.08957
# 10 chr10:119691120:G:C 119691120 C ADD 1644 1.426 0.1185 1.13 1.798 2.992 0.002776 G C G 0.354873
# 2 chr2:178395655:A:G 178395655 A ADD 1643 0.7096 0.1157 0.5657 0.8902 -2.966 0.003017 A G G -0.343054
# 10 chr10:119610713:G:A 119610713 A ADD 1645 2.82 0.3518 1.415 5.619 2.947 0.003205 G A G 1.03674
# 10 chr10:119877390:G:A 119877390 A ADD 1645 2.966 0.3709 1.434 6.136 2.932 0.003371 G A G 1.08721


# CCSS EXP
awk '{split($2, a, ":"); print $0, a[3], a[4]}' ccss_exp_results.assoc.logistic.clean.Psorted \
| awk '{if($4==$13) A2=$14; else A2=$13; BETA=log($7); print $0, A2, BETA}' \
| awk 'BEGIN{print "CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA"}FNR>1{print}' \
> ccss_exp_results.assoc.logistic.clean.Psorted.formetal
# CCSS ORG
awk '{split($2, a, ":"); print $0, a[3], a[4]}' ccss_org_results.assoc.logistic.clean.Psorted \
| awk '{if($4==$13) A2=$14; else A2=$13; BETA=log($7); print $0, A2, BETA}' \
| awk 'BEGIN{print "CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA"}FNR>1{print}' \
> ccss_org_results.assoc.logistic.clean.Psorted.formetal

## Identify variants that match between the 3 datasets, based on chr:pos and then based on both A1 and A2
awk 'NR==FNR{a[$1":"$3]=$0;next}($1":"$3 in a){print $0, a[$1":"$3]}' ccss_org_results.assoc.logistic.clean.Psorted.formetal \
sjlife_results.assoc.logistic.clean.Psorted.formetal \
| awk '($4==$20 || $4==$31) && ($15==$20 || $15==$31)' \
|  awk 'NR==FNR{a[$1":"$3]=$0;next}($1":"$3 in a){print $0, a[$1":"$3]}' - ccss_exp_results.assoc.logistic.clean.Psorted.formetal \
| awk '($4==$36 || $4==$47) && ($15==$36 || $15==$47)' \
> common_variants_among_3_datasets.txt
# There were 3 variants that were duplicates and because they were all not-significant, retain only unique records
# 10 chr10:119528844:TAATA:T 119528844 TAATA ADD 1452 1.306 0.2085 0.8677 1.965 1.279 0.2008 TAATA T T 0.266969 10 chr10:119528844:T:TAATA 119528844 TAATA ADD 1636 1.002 0.2146 0.6582 1.527 0.0115 0.9908 T TAATA T 0.001998 10 chr10:119528844:T:TAATA 119528844 TAATA ADD 2599 0.8836 0.2593 0.5316 1.469 -0.4773 0.6331 T TAATA T -0.123751
# 10 chr10:119875979:TA:T 119875979 T ADD 1374 1.374 0.3574 0.6821 2.769 0.8896 0.3737 TA T TA 0.317726 10 chr10:119875979:TA:T 119875979 T ADD 1501 0.877 0.2612 0.5256 1.463 -0.5024 0.6154 TA T TA -0.131248 10 chr10:119875979:TA:T 119875979 T ADD 3126 1.038 0.2266 0.666 1.619 0.1659 0.8682 TA T TA 0.0372958
# 10 chr10:119614280:AAAAT:A 119614280 A ADD 1456 0.8726 0.4685 0.3484 2.186 -0.2908 0.7712 AAAAT A AAAAT -0.136278 10 chr10:119614280:A:AAAAT;chr10:119614280:A:AAAATAAATAAAT 119614280 AAAAT ADD 1635 0.9481 0.1573 0.6965 1.291 -0.3389 0.7347 A AAAAT;chr10 A -0.0532953 10 chr10:119614280:AAAAT:A 119614280 A ADD 2976 1.102 0.4915 0.4207 2.888 0.198 0.8431 AAAAT A AAAAT 0.0971267
# Then prepare final files for each dataset for METAL
awk '!a[$1":"$3]++' common_variants_among_3_datasets.txt | cut -d' ' -f1-16 | awk '{if(FNR>1) $2=$1":"$3; print}' \
> ccss_exp_results.assoc.logistic.clean.Psorted.formetal.common
awk '!a[$1":"$3]++' common_variants_among_3_datasets.txt | cut -d' ' -f17-32 | awk '{if(FNR>1) $2=$1":"$3; print}' \
> sjlife_results.assoc.logistic.clean.Psorted.formetal.common
awk '!a[$1":"$3]++' common_variants_among_3_datasets.txt | cut -d' ' -f33-48 | awk '{if(FNR>1) $2=$1":"$3; print}' \
> ccss_org_results.assoc.logistic.clean.Psorted.formetal.common

## Prepare the config file for METAL
# run_metal.script

## Run METAL
~/bin/generic-metal/metal run_metal.script

## A total of 2601 variants present in all 3 datasets were meta-analyzed; process the results further
