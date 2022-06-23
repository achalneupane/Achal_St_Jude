#!/bin/bash

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs

mkdir -p prs_out

study=$1

# Subset PRS data for each study
awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt > prs_out/ALL_Cancers_PRS_data.txt_$study
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/ALL_Cancers_PRS_data.txt_$study | wc -l
# Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/ALL_Cancers_PRS_data.txt_$study plink_data/sjlife_all_PRS_all.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match
# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/ALL_Cancers_PRS_data.txt_$study plink_data/sjlife_all_PRS_all.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match \
> prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

wc -l prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

grep -vw chr1:113903258:G:T prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final
# Check for duplicate variants
# awk 'a[$1":"$3]++' prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_duplicates
# Drop one from the duplicate; check for the allele frequency first, and get rid of the rare variant keeping the common one
# egrep -vw 'chr9:108126198:G:A|chr3:30641447:G:C' prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_uniq
# grep -vw chr9:108126198:G:A prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_uniq
# mv prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_uniq prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match
# Harmonize no direct match alleles
module load R
Rscript harmonize_alleles.R prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final
# Update the alleles
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# Extract study-specific variants
awk '{print $2}' prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match_to_extract.txt
module load plink/1.90b
plink --bfile plink_data/sjlife_all_PRS_all --extract prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
plink --bfile plink_data/sjlife_all_PRS_all --extract prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized
plink --bfile prs_out/${study}_direct_match --bmerge prs_out/${study}_harmonized --make-bed --out prs_out/$study
# Update variant names
awk '{print $2, $1":"$4}' prs_out/${study}.bim > prs_out/${study}_update_variantnames
plink --bfile prs_out/$study --update-name prs_out/${study}_update_variantnames --make-bed --out prs_out/${study}_varname_updated
# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out/ALL_Cancers_PRS_data.txt_$study > prs_out/${study}.prsweight
# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs


# ALL_Vijayakrishnan	15
# Allman_African_Breast	75	74
# Allman_Hispanic_Breast	71 69
# Khera_2018_Breast	77	76
# Mavaddat_2015_ER_NEG_Breast	77	76
# Mavaddat_2015_ER_OVERALL_Breast	77
# Mavaddat_2015_ER_POS_Breast	77	76
# Mavaddat_2019_ER_NEG_Breast	313
# Mavaddat_2019_ER_OVERALL_Breast	313
# Mavaddat_2019_ER_POS_Breast	313
# Meningioma_Claus	1
# Meningioma_Dobbins	2
# MichiganWeb_ER_NEG_Breast	79
# MichiganWeb_ER_OVERALL_Breast	1120348
# MichiganWeb_ER_POS_Breast	1119079
# Pleiotropy_Bi_directional	21
# Pleiotropy_Meta_analysis	21
# Pleiotropy_One_cohort	9
# Pleiotropy_One_directional	137
# Pleiotropy_PRSWEB	179
# Pleiotropy_Replication_prior_studies	308
# Sarcoma_Machiela	6
# Wang_African_Breast	98



# # harmonize_allele.R
# # Process the variants to match their alleles

# # Read the input file including variants with no direct match of their alleles
# args = commandArgs(trailingOnly = TRUE)
# dat = read.table(args[1], header = FALSE, stringsAsFactors = FALSE)
# # dat = read.table("Z:/ResearchHome/ClusterHome/ysapkota/Work/CAD_PRS/y", header = FALSE, stringsAsFactors = FALSE)

# ## Function to flip alleles
# flip_alleles = function(x){
#   if(x=="A"){
#     y="T"
#   } else if (x=="T"){
#     y="A"
#   } else if (x=="C"){
#     y="G"
#   } else if (x=="G"){
#     y="C"
#   }
#   return(y)
# }

# # Process each variant to check the alleles
# dat.out = NULL
# for (i in 1:nrow(dat)){
#   chr=dat$V1[i]
#   pos = dat$V3[i]
#   variant = dat$V2[i]
#   khera_a1 = dat$V6[i]
#   khera_a2 = dat$V7[i]
#   # khera_weight = dat$V9[i]
#   wgs_a1 = dat$V4[i]
#   wgs_a2 = dat$V5[i]
#   # First find out if the wgs_a1 have more than one character
#   if(nchar(wgs_a1)>1){
#     wgs_a1_first = substr(wgs_a1, 1, 1)
#     wgs_a1_last = substr(wgs_a1, nchar(wgs_a1), nchar(wgs_a1))
#     wgs_a1_changed = ifelse(wgs_a1_first==wgs_a2, wgs_a1_last, wgs_a1_first)
#   } else {
#     wgs_a1_changed = wgs_a1
#   }
#   # Then do the same for wgs_a2 allele
#   if (nchar(wgs_a2)>1){
#     wgs_a2_first = substr(wgs_a2, 1, 1)
#     wgs_a2_last = substr(wgs_a2, nchar(wgs_a2), nchar(wgs_a2))
#     wgs_a2_changed = ifelse(wgs_a2_first==wgs_a1, wgs_a2_last, wgs_a2_first)
#   } else {
#     wgs_a2_changed = wgs_a2
#   }
#   # Now check if the changed wgs alleles match with those from Khera et al
#   # No alleles flipped
#   if ((wgs_a1_changed == khera_a1 & wgs_a2_changed == khera_a2) | (wgs_a1_changed == khera_a2 & wgs_a2_changed == khera_a1)){
#     wgs_a1_new = wgs_a1_changed; wgs_a2_new = wgs_a2_changed; match=1
#   } else if((flip_alleles(wgs_a1_changed) == khera_a1 & wgs_a2_changed == khera_a2) | (flip_alleles(wgs_a1_changed) == khera_a2 & wgs_a2_changed == khera_a1)) { # only a1 flipped
#     wgs_a1_new = flip_alleles(wgs_a1_changed); wgs_a2_new = wgs_a2_changed; match=1
#   } else if ((wgs_a1_changed == khera_a1 & flip_alleles(wgs_a2_changed) == khera_a2) | (wgs_a1_changed == khera_a2 & flip_alleles(wgs_a2_changed) == khera_a1)) { # only a2 flipped
#     wgs_a1_new = wgs_a1_changed; wgs_a2_new = flip_alleles(wgs_a2_changed); match=1
#   } else if ((flip_alleles(wgs_a1_changed) == khera_a1 & flip_alleles(wgs_a2_changed) == khera_a2) | (flip_alleles(wgs_a1_changed) == khera_a2 & flip_alleles(wgs_a2_changed) == khera_a1)){ # both a1 and a2 flipped
#     wgs_a1_new = flip_alleles(wgs_a1_changed); wgs_a2_new = flip_alleles(wgs_a2_changed); match=1
#   } else {
#     wgs_a1_new = wgs_a1_changed; wgs_a2_new = wgs_a2_changed; match=0
#   }
#   dat.out = rbind(dat.out, data.frame(chr, pos, variant, khera_a1, khera_a2, wgs_a1, wgs_a2, wgs_a1_new, wgs_a2_new, match))
# }

# # Write data to disc
# write.table(dat.out, paste0(args[1], "_alleles_harmonized"), row.names = FALSE, quote = FALSE)