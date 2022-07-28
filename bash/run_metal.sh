#!/bin/bash

## Process PLINK result files for meta-analysis using METAL; add A2 and BETA
# SJLIFE
awk '{split($2, a, ":"); print $0, a[3], a[4]}' sjlife_results.assoc.logistic.clean.Psorted \
| awk '{if($4==$13) A2=$14; else A2=$13; BETA=log($7); print $0, A2, BETA}' \
| awk 'BEGIN{print "CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA"}FNR>1{print}' \
> sjlife_results.assoc.logistic.clean.Psorted.formetal
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
