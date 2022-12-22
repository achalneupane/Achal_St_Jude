#!/bin/bash

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs

mkdir -p prs_out

study=$1
study="Basal_cell_carcinoma_PRSWeb"
study="Squamous_cell_carcinoma_PRSWeb"
study="Mavaddat_2019_ER_NEG_Breast"
# Subset PRS data for each study
## remove chr1:145902073|chr4:57426897|chr6:114515866 from MichiganWeb_ER_OVERALL_Breast
## remove chr4:57426897|chr6:114515866 from MichiganWeb_ER_POS_Breast
awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt > prs_out/ALL_Cancers_PRS_data.txt_${study}
# awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt | egrep -v '145902073|57426897|114515866|129989587' > prs_out/ALL_Cancers_PRS_data.txt_${study} # MichiganWeb_ER_OVERALL_Breast
# awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt | egrep -v '57426897|114515866' > prs_out/ALL_Cancers_PRS_data.txt_${study} # MichiganWeb_ER_POS_Breast
# awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt | grep -v 30641447 > prs_out/ALL_Cancers_PRS_data.txt_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/ALL_Cancers_PRS_data.txt_$study | wc -l
# Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/ALL_Cancers_PRS_data.txt_$study plink_data/sjlife_all_PRS_all_final_v3.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match
# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/ALL_Cancers_PRS_data.txt_$study plink_data/sjlife_all_PRS_all_final_v3.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match \
> prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

wc -l prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

# grep -vw chr1:113903258:G:T prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final
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
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final_alleles_harmonized > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final_alleles_harmonized | grep -vw chr9:108126198:G:A > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# grep -v chr1:145902073:G:GA prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt > prs_out/t1
# mv prs_out/t1 prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# Extract study-specific variants
awk '{print $2}' prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match_to_extract.txt
module load plink/1.90b
plink --bfile plink_data/sjlife_all_PRS_all_final_v3 --extract prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
plink --bfile plink_data/sjlife_all_PRS_all_final_v3 --extract prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized
plink --bfile prs_out/${study}_direct_match --bmerge prs_out/${study}_harmonized --make-bed --out prs_out/$study
# Update variant names
awk '{print $2, $1":"$4}' prs_out/${study}.bim > prs_out/${study}_update_variantnames
plink --bfile prs_out/$study --update-name prs_out/${study}_update_variantnames --make-bed --out prs_out/${study}_varname_updated
# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out/ALL_Cancers_PRS_data.txt_$study > prs_out/${study}.prsweight
# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs

## _Significant are the ones with: grep -v N
# ALL_Vijayakrishnan	15 15 # OK
# Allman_African_Breast	75	75 # OK; dropped chr9:108126198:G:A duplicate
# Allman_Hispanic_Breast	71 70 # chr3:30641447:G:G in PRS file; removed
# Khera_2018_Breast	77	77 # OK
# Mavaddat_2015_ER_NEG_Breast	77	77 # OK
# Mavaddat_2015_ER_OVERALL_Breast	77 77 # OK
# Mavaddat_2015_ER_POS_Breast	77	77 # OK
# Mavaddat_2019_ER_NEG_Breast	313 306
# Mavaddat_2019_ER_OVERALL_Breast	313 306
# Mavaddat_2019_ER_POS_Breast	313 306
# Meningioma_Claus	1 1 #OK
# Meningioma_Dobbins	2 2 # OK
# Meningioma 3 3 # OK
# MichiganWeb_ER_NEG_Breast	79 79
# MichiganWeb_ER_OVERALL_Breast	1120348 1119287
# MichiganWeb_ER_POS_Breast	1119079 1118049
# Pleiotropy_Bi_directional_Increasing	21 21 # OK
# Pleiotropy_Bi_directional_Increasing_Significant	15 15 # OK
# Pleiotropy_Bi_directional_Decreasing	21 21 # OK
# Pleiotropy_Bi_directional_Decreasing_Significant 15 15 # OK
# Pleiotropy_Meta_analysis	21 20 # OK; removed one duplicate variant
# Pleiotropy_One_cohort	9 9 # OK
# Pleiotropy_One_directional	137 137 # OK
# Pleiotropy_One_directional_Significant	85 85 # OK
# Pleiotropy_PRSWEB	179 179 # OK
# Pleiotropy_Replication_prior_studies	308 308 # OK
# Sarcoma_Machiela	6 6 # OK
# Wang_African_Breast	98 98
# sjlife_thyroid.profile 12 12




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






# # MichiganWeb_ER_OVERALL_Breast; 
# Extract duplicate vars
plink --bfile prs_out/MichiganWeb_ER_OVERALL_Breast_varname_updated --list-duplicate-vars ids-only suppress-first
# Error: Duplicate ID : 1:145902073|4:57426897; so removing them with maf of lower values from prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt 


###################################
## Find those that did not match ##
###################################
## Additional notes:
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Allman_African_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_Allman_African_Breast
Allman_African_Breast
12      114691  G       A       0.10436 Allman_African_Breast   Breast  Y
# chr12_GL383549v1_alt	114691	# replaced with the correct position: chr12:28002147

awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Allman_Hispanic_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_Allman_Hispanic_Breast
Allman_Hispanic_Breast
3       30641447        G       G       0.06765865      Allman_Hispanic_Breast  Breast  Y
12      114691  C       T       0.3435897       Allman_Hispanic_Breast  Breast  Y
# chr12_GL383549v1_alt	114691	# replaced with the correct position: chr12:28002147
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Khera_2018_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_Khera_2018_Breast
Khera_2018_Breast
12      114691  G       A       0.147456469     Khera_2018_Breast       Breast  Y
# chr12_GL383549v1_alt	114691	# replaced with the correct position: chr12:28002147
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Mavaddat_2015_ER_NEG_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_Mavaddat_2015_ER_NEG_Breast
Mavaddat_2015_ER_NEG_Breast
12      114691  G       A       0.176617854     Mavaddat_2015_ER_NEG_Breast     Breast  Y
# chr12_GL383549v1_alt	114691	# replaced with the correct position: chr12:28002147
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Mavaddat_2015_ER_POS_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_Mavaddat_2015_ER_POS_Breast
Mavaddat_2015_ER_POS_Breast
12      114691  G       A       0.126243727     Mavaddat_2015_ER_POS_Breast     Breast  Y
# chr12_GL383549v1_alt	114691	# replaced with the correct position: chr12:28002147
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Mavaddat_2019_ER_NEG_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_Mavaddat_2019_ER_NEG_Breast
Mavaddat_2019_ER_NEG_Breast
# 1       145830809       CT      C       0.0126  Mavaddat_2019_ER_NEG_Breast     Breast  Y # replaced with position 145830798; no matching allele
# 1       204533386       TTCTGAAACAGGG   T       0.1345  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 2       217091173       G       GA      0.0558  Mavaddat_2019_ER_NEG_Breast     Breast  Y ## No matching allele
# 3       55936749        AT      A       0.0586  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 3       63901773        T       TTG     0.043   Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 4       125831837       AAT     A       0.0638  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 4       83448971        TA      TAA     0.0489  Mavaddat_2019_ER_NEG_Breast     Breast  Y ## No matching allele
# 5       1296140 AG      A       0.1056  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 5       58945885        T       C       0.0408  Mavaddat_2019_ER_NEG_Breast     Breast  Y ## No match
# 5       79885172        G       GA      0.0804  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 6       130020583       C       CT      0.0804  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 6       87094101        T       C       0.0678  Mavaddat_2019_ER_NEG_Breast     Breast  Y ## Not found; confirmed
# 7       91829875        A       ATT     0.0486  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 8       127201316       CA      C       0.04    Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 11      108396675       CA      C       0.0629  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 12      82670416        G       GA      0.0717  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 17      30841059        T       G       0.0604  Mavaddat_2019_ER_NEG_Breast     Breast  Y ## Not found; confirmed, but have variant for 30841058 that do not match alleles
# 19      19406245        C       CGGGCG  0.0577  Mavaddat_2019_ER_NEG_Breast     Breast  Y
# 22      38187308        AAAAG   AAAAGAAAG       0.0079  Mavaddat_2019_ER_NEG_Breast     Breast  Y # No matching allele
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' MichiganWeb_ER_NEG_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_MichiganWeb_ER_NEG_Breast
MichiganWeb_ER_NEG_Breast
# 1       145790095       T       C       0.0512933143875503      MichiganWeb_ER_NEG_Breast       Breast  Y ## Replaced with chr1:145790097
Pleiotropy_Meta_analysis has two duplicate variants from different studies. We decided to use the one with the beta -0.916290732
# chr17	7668434	G	T	-0.916290732
# chr17	7668434	G	T	-0.400477567
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Pleiotropy_One_directional_varname_updated.bim  ALL_Cancers_PRS_data.txt_Pleiotropy_One_directional
Pleiotropy_One_directional
# 17      45720982        T       TTTG    0.075801713     Pleiotropy_One_directional      Pleiotropy      Y
# 10      46037695        T       C       0.104250021     Pleiotropy_One_directional      Pleiotropy      Y ## replaced with chr 10:46037697
# 2       200810808       C       CG      0.06720875      Pleiotropy_One_directional      Pleiotropy      N
# 22      38185855        T       TC      0.048140375     Pleiotropy_One_directional      Pleiotropy      N 
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Pleiotropy_Replication_prior_studies_varname_updated.bim  ALL_Cancers_PRS_data.txt_Pleiotropy_Replication_prior_studies
Pleiotropy_Replication_prior_studies
# 10      46046324        C       T       0.222743471     Pleiotropy_Replication_prior_studies    Pleiotropy      Y ## replaced with chr10:46046326
# 9       108131271       C       CT      0.095630231     Pleiotropy_Replication_prior_studies    Pleiotropy      Y
# 3       128027236       AT      A       0.081579987     Pleiotropy_Replication_prior_studies    Pleiotropy      Y 
# 10      46037695        T       C       0.099709844     Pleiotropy_Replication_prior_studies    Pleiotropy      Y ## replaced with chr 10:46037697
# 14      68174379        CA      C       0.080657903     Pleiotropy_Replication_prior_studies    Pleiotropy      Y
# 2       241237728       CAT     C       0.083605584     Pleiotropy_Replication_prior_studies    Pleiotropy      Y
# 1       150549995       C       CTG     0.079043207     Pleiotropy_Replication_prior_studies    Pleiotropy      Y
# 6       21159817        T       TA      0.105249411     Pleiotropy_Replication_prior_studies    Pleiotropy      Y
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' Wang_African_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_Wang_African_Breast
Wang_African_Breast
# 1       145790095       T       C       0.040821995     Wang_African_Breast     Breast  Y ## Replaced with chr1:145790097

awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' MichiganWeb_ER_OVERALL_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_MichiganWeb_ER_OVERALL_Breast
MichiganWeb_ER_OVERALL_Breast
## These are the variants that need to be re-checked
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' MichiganWeb_ER_OVERALL_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_MichiganWeb_ER_OVERALL_Breast | \
awk '{print "chr"$1"\t"$2}' > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/MichiganWeb_ER_OVERALL_Breast_reconfirm_variants.txt
## Extract those that needs to be rechecked
awk -v OFS="\t" '{sub(/_.*/,"",$1)}1' Michigan_OVERALL_GrCh38.bed > Michigan_OVERALL_GrCh38_chr_only.bed
awk 'NR==FNR{a[$1,$2];next} ($1,$2) in a' MichiganWeb_ER_OVERALL_Breast_reconfirm_variants.txt Michigan_OVERALL_GrCh38_chr_only.bed > MichiganWeb_ER_OVERALL_Breast_reconfirm_variants_GR37_Gr38.txt
# Bed file
awk 'BEGIN {FS=OFS="\t"}
{
  split($4,a,":")
  print $0 OFS a[1] ":" a[2]-1 "-" a[2]
}
' MichiganWeb_ER_OVERALL_Breast_reconfirm_variants_GR37_Gr38.txt > MichiganWeb_ER_OVERALL_Breast_reconfirm_variants_GR37_Gr38_for_bed.txt

awk '{print $4 $5}' MichiganWeb_ER_OVERALL_Breast_reconfirm_variants_GR37_Gr38_for_bed.txt

awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' MichiganWeb_ER_POS_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_MichiganWeb_ER_POS_Breast
MichiganWeb_ER_POS_Breast
## These are the variants that need to be re-checked
awk 'NR==FNR{a[$2]; next} !(($1":"$2) in a)' MichiganWeb_ER_POS_Breast_varname_updated.bim  ALL_Cancers_PRS_data.txt_MichiganWeb_ER_POS_Breast | \
awk '{print "chr"$1"\t"$2}'> /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/MichiganWeb_ER_POS_Breast_reconfirm_variants.txt
## Extract those that needs to be rechecked
awk -v OFS="\t" '{sub(/_.*/,"",$1)}1' Michigan_ER_POS_GrCh38.bed > Michigan_ER_POS_GrCh38_chr_only.bed
awk 'NR==FNR{a[$1,$2];next} ($1,$2) in a' MichiganWeb_ER_POS_Breast_reconfirm_variants.txt Michigan_ER_POS_GrCh38_chr_only.bed > MichiganWeb_ER_POS_Breast_reconfirm_variants_GR37_Gr38.txt

# Bed file
awk 'BEGIN {FS=OFS="\t"}
{
  split($4,a,":")
  print $0 OFS a[1] ":" a[2]-1 "-" a[2]
}
' MichiganWeb_ER_POS_Breast_reconfirm_variants_GR37_Gr38.txt > MichiganWeb_ER_POS_Breast_reconfirm_variants_GR37_Gr38_for_bed.txt



# ## Re-doing the lift-over for the 1057 varinats from Michigan_ER_OVERALL
# /home/aneupane/liftover/liftOver Michigan_ER_OVERALL_liftover_reconfirm.bed /home/aneupane/liftover/hg19ToHg38.over.chain Michigan_OVERALL_GrCh38_1057_vars_reconfirm.bed Michigan_OVERALL_GrCh38_1057_vars_reconfirm_unmapped.bed



MichiganWeb_ER_OVERALL_Breast
[aneupane@noderome157 prs]$ grep -w 57426897  prs_out/${study}.bim
4       chr4:57426897:T:TG      0       57426897        G       T
4       chr4:57426897:T:TGTG    0       57426897        G       T
4       chr4:57426897:T:TGTGTG  0       57426897        G       T
4       chr4:57426897:T:TGTGTGTG        0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTG      0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTGTG    0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTGTGTG  0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTGTGTGTG        0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTGTGTGTGTGTG    0       57426897        G       T
4       chr4:57426897:T:TTG     0       57426897        G       T
4       chr4:57426897:T:TTGTG   0       57426897        G       T
4       chr4:57426897:T:TTGTGTG 0       57426897        G       T
4       chr4:57426897:T:TTGTGTGTG       0       57426897        G       T
4       chr4:57426897:TTG:T     0       57426897        T       G
4       chr4:57426897:TTGTG:T   0       57426897        T       G
4       chr4:57426897:TTGTGTG:T 0       57426897        T       G
4       chr4:57426897:TTGTGTGTG:T       0       57426897        T       G

[aneupane@noderome157 prs]$ grep -w 57426897  plink_data/plink.frq

   4    chr4:57426897:T:TG   TG    T    0.0006664     9004
   4    chr4:57426897:T:TGTG TGTG    T     0.004331     9004
   4    chr4:57426897:T:TGTGTG TGTGTG    T     0.009662     9004
   4    chr4:57426897:T:TGTGTGTG TGTGTGTG    T       0.1308     9004
   4    chr4:57426897:T:TGTGTGTGTG TGTGTGTGTG    T       0.1825     9004
   4    chr4:57426897:T:TGTGTGTGTGTG TGTGTGTGTGTG    T       0.0743     9004
   4    chr4:57426897:T:TGTGTGTGTGTGTG TGTGTGTGTGTGTG    T     0.004442     9004
   4    chr4:57426897:T:TGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTG    T    0.0005553     9004
   4    chr4:57426897:T:TGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTG    T    0.0002221     9004
   4    chr4:57426897:T:TTG  TTG    T       0.1388     9004
   4    chr4:57426897:T:TTGTG TTGTG    T      0.01166     9004
   4    chr4:57426897:T:TTGTGTG TTGTGTG    T     0.006219     9004
   4    chr4:57426897:T:TTGTGTGTG TTGTGTGTG    T     0.001111     9004
   4    chr4:57426897:TTG:T    T  TTG      0.06319     9004
   4    chr4:57426897:TTGTG:T    T TTGTG      0.09429     9004
   4    chr4:57426897:TTGTGTG:T    T TTGTGTG     0.005886     9004
   4    chr4:57426897:TTGTGTGTG:<*:DEL> <*:DEL> TTGTGTGTG    0.0003332     9004
   4    chr4:57426897:TTGTGTGTG:T    T TTGTGTGTG    0.0003332     9004

# # Mavaddat_2019_ER_*_Breast: 4 no match; 15 not in WGS= 19 in total variant ID mismatch
# chr pos variant khera_a1 khera_a2 wgs_a1 wgs_a2 wgs_a1_new wgs_a2_new match
# 4 83448971 chr4:83448971:TA:T TA TAA TA T A T 0
# 7 91829875 chr7:91829875:AT:A A ATT A AT A T 0
# 22 38187308 chr22:38187308:AAAAG:A AAAAG AAAAGAAAG A AAAAG A G 0
# 22 38187308 chr22:38187308:AAAAGAAAG:A AAAAG AAAAGAAAG A AAAAGAAAG A G 0




4       chr4:57426897:T:TG      0       57426897        G       T
4       chr4:57426897:T:TGTG    0       57426897        G       T
4       chr4:57426897:T:TGTGTG  0       57426897        G       T
4       chr4:57426897:T:TGTGTGTG        0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTG      0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTGTG    0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTGTGTG  0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTGTGTGTG        0       57426897        G       T
4       chr4:57426897:T:TGTGTGTGTGTGTGTGTGTG    0       57426897        G       T
4       chr4:57426897:T:TTG     0       57426897        G       T
4       chr4:57426897:T:TTGTG   0       57426897        G       T
4       chr4:57426897:T:TTGTGTG 0       57426897        G       T
4       chr4:57426897:T:TTGTGTGTG       0       57426897        G       T
4       chr4:57426897:TTG:T     0       57426897        T       G
4       chr4:57426897:TTGTG:T   0       57426897        T       G
4       chr4:57426897:TTGTGTG:T 0       57426897        T       G
4       chr4:57426897:TTGTGTGTG:T       0       57426897        T       G

Aleles_in_dbSNP    sjlife_EUR_maf    1kg_EUR_maf
rs1851999
chr4:57426897:T:A
chr4:57426897:T:G    				0.4
rs1553882276
chr4:57426897:T:TG    0.0001455
chr4:57426897:T:TGTG    0.00422
chr4:57426897:T:TGTGTG    0.00975
chr4:57426897:T:TGTGTGTG    0.1342
chr4:57426897:T:TGTGTGTGTG    0.2125
chr4:57426897:T:TGTGTGTGTGTG    0.08411
chr4:57426897:T:TGTGTGTGTGTGTG    0.00422
chr4:57426897:T:TGTGTGTGTGTGTGTG    0.0005821
rs1717025271
chr4:57426897:T:TTGTGTA

# extract from SJLIFE
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/plink_data
KEEP=/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_EUR_Per_PCA.txt
plink --bfile sjlife_all_PRS_all_final --keep ${KEEP} --make-bed --freq --out test_freq_Eur

# extract from 1000 genome
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/1kGP
plink --bfile 1000genomes_merged --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/1KG_EUR.txt --make-bed --out test_freq_Eur
awk '{print "chr"$1"\tchr"$1":"$4":"$5":"$6"\t"$3"\t"$4"\t"$6"\t"$5}' test_freq_Eur.bim> test_freq_Eur.bim2

plink --bfile test_freq_Eur --freq --out test_freq_Eur_freq


##################################################
## Part 2: CCSS org overlapping PRS from SJLIFE ##
##################################################
## Since CCSS_original was missing 72 variants in the GWAS data, I am re-calculating PRS for breast cancer with just 241 variants; Meningioma with just 2 overlapping variants; Sarcoma with 5 overlapping variants; pleiotropy with 172 overlapping variants.
## Use R Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/CCSS/ccss_org_bed.R (Part 2) to extract variants
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs

study="Meningioma_from_variants_also_in_CCSS_org"
awk -v study=$study '$6==study' ALL_Cancers_PRS_data_in_CCSS_org.txt > prs_out/ALL_Cancers_PRS_data.txt_${study}
# ## Meningioma
## grep "Meningioma" ALL_Cancers_PRS_data_in_CCSS_org.txt > prs_out/ALL_Cancers_PRS_data.txt_${study}

# awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt | egrep -v '145902073|57426897|114515866|129989587' > prs_out/ALL_Cancers_PRS_data.txt_${study} # MichiganWeb_ER_OVERALL_Breast
# awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt | egrep -v '57426897|114515866' > prs_out/ALL_Cancers_PRS_data.txt_${study} # MichiganWeb_ER_POS_Breast
# awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt | grep -v 30641447 > prs_out/ALL_Cancers_PRS_data.txt_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/ALL_Cancers_PRS_data.txt_$study | wc -l
# Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/ALL_Cancers_PRS_data.txt_$study plink_data/sjlife_all_PRS_all_final_v3.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match
# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/ALL_Cancers_PRS_data.txt_$study plink_data/sjlife_all_PRS_all_final_v3.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match \
> prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

wc -l prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

# check how many in PRS and direct match
wc -l prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match
grep $study ALL_Cancers_PRS_data_in_CCSS_org.txt| wc -l

# grep -vw chr1:113903258:G:T prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final
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
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final_alleles_harmonized > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final_alleles_harmonized | grep -vw chr9:108126198:G:A > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# grep -v chr1:145902073:G:GA prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt > prs_out/t1
# mv prs_out/t1 prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# Extract study-specific variants
awk '{print $2}' prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match_to_extract.txt
module load plink/1.90b
plink --bfile plink_data/sjlife_all_PRS_all_final_v3 --extract prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
plink --bfile plink_data/sjlife_all_PRS_all_final_v3 --extract prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized
plink --bfile prs_out/${study}_direct_match --bmerge prs_out/${study}_harmonized --make-bed --out prs_out/$study
# Update variant names
awk '{print $2, $1":"$4}' prs_out/${study}.bim > prs_out/${study}_update_variantnames
plink --bfile prs_out/$study --update-name prs_out/${study}_update_variantnames --make-bed --out prs_out/${study}_varname_updated
# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out/ALL_Cancers_PRS_data.txt_$study > prs_out/${study}.prsweight
# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs





CHR  POS REF EA BETA
1 33333 G A 0.03



1       chr1:33333        0       33333 A       G



