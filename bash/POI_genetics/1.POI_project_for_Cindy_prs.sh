# Hi Achal,

# Below are three studies with Cindy which will require P/LP on different sets of genes. For POI study, we need to calculate a PRS based on published study. Can you please include in your list of tasks to do and work on them once we address the issues with P/LP variants? You can calculate the PRS when you get time, however.

# Thanks,

# Yadav

# 2 - POI genetics
# 2a - Achal(?) will compute PRS from GWS SNPs from Ruth et al and P/LPs from Ke et al and Shekari et al (would end of Sep. be okay for this?); analyst Kenneth will work on phenotype cleaning and review
# 2b - POI GWAS in survivors (GRA Aparna will work on this)
# 2c - P/LP-PRS analysis (analyst Kenneth will work on this)

# 3 - BCC genetics
# 3a - Achal will provide updated P/LP for cancer and BCC as appropriate (e.g., SJCPG60, Kim172, gene panels), please provide by end of Sep(?)
# 3b - Clinical risk prediction + BCC genetics update (analyst Kenneth will finish this)
# 3c - BCC GWAS in survivors (Cindy will work on this)
 
# 4 - ACMG
# 4a - Achal will provide P/LP by gene, condition grouping, etc. (done by end of Oct?)
# 4b - Preliminary descriptive summary: analyst Kenneth can complete this because he will be doing something similar for 3b (Yadav, do you want Achal or someone else to do this?)
# 4c - Other analysis: Kendrick will coordinate


## Extract variants for PRS
# load modules
module load bcftools
module load plink/1.90b


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/plink_data/
################################################
## create plink files for sjlife and CCSS exp ##
################################################
module load plink/1.90b
for CHR in {1..22} X Y; do
echo "Doing Chr ${CHR}"
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/preQC/Survivor_WGS.GATK4180.hg38_renamed_chr${CHR}.decomposed.preqc --extract ../wantedSNP_GRCh38.txt --keep-allele-order --make-bed --out Ruth_etal_PRS_chr${CHR}
done


ls Ruth_etal_PRS_chr*.bed | sed 's/\..*//' | sort -u > merge_list.txt
plink --merge-list merge_list.txt -keep-allele-order --make-bed --out Ruth_etal_PRS_merged_Survivor_WGS


## After running R code PRS_scores.R; extract variants
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/
ln -s  ../POI_meta_GRCh38.dat .
ln -s  ../POI_meta_GRCh37.dat .


ln -s ../plink_data/Ruth_etal_PRS_merged_Survivor_WGS.* .

# Append chr in the first column
awk '{print "chr"$1, $2, $3, $4, $5, $6}' Ruth_etal_PRS_merged_Survivor_WGS.bim > tmp_bim
mv tmp_bim Ruth_etal_PRS_merged_Survivor_WGS.bim

## Note: SJLIFE and CCSS_exp have all variants; CCSS_org is missing 9 variants in total
## PRS; studies: PGS000356 (179 vars) PGS000454 (27 vars) PGS003416 (462 vars) ST6 (110 vars)
mkdir -p prs_out_wgs_survivor
study="POI_META"
# Subset PRS data for each study
awk -v study=$study '$6==study' POI_meta_GRCh38.dat > prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study} | wc -l
## Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study} Ruth_etal_PRS_merged_Survivor_WGS.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_direct_match

# (base) [aneupane@splprhpc12 prs]$ wc -l prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_direct_match
# 283 prs_out_wgs_survivor/POI_meta_GRCh38.dat_POI_META_direct_match

# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study} Ruth_etal_PRS_merged_Survivor_WGS.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v '*' > prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_direct_match prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match \
> prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match_final

wc -l prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match_final
# 6 prs_out_wgs_survivor/POI_meta_GRCh38.dat_POI_META_no_direct_match_final

## Harmonize alleles
module load R/4.2.2-rhel8
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match_final

## Update the alleles; It checks if the last field ($NF) in each line is equal to 1; if true, then { print $3, $6, $7, $8, $9 } so we can update the plink file for variants with no direct match
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match_final_alleles_harmonized > prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt


# Extract study-specific variants
awk '{print $2}' prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_direct_match > prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_direct_match_to_extract.txt
module load plink/1.90b

## If there are no indirect match, just make this commented first line: --out prs_out_wgs_survivor/${study} and skip lines with **
plink --bfile Ruth_etal_PRS_merged_Survivor_WGS --extract prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out_wgs_survivor/${study}_direct_match
# plink --bfile Ruth_etal_PRS_merged_Survivor_WGS --extract prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out_wgs_survivor/${study}
plink --bfile Ruth_etal_PRS_merged_Survivor_WGS --extract prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out_wgs_survivor/${study}_harmonized_tmp ##**
# This command above extracts and updates alleles from 
# (base) [aneupane@splprhpc12 prs]$ cat prs_out_wgs_survivor/POI_META_harmonized.bim
# 13      chr13:40690558:A:AT     0       40690558        T       A
# 13      chr13:40690558:A:ATTTT  0       40690558        T       A
# 13      chr13:40690558:A:ATTTTT 0       40690558        T       A
# 13      chr13:40690558:A:ATTTTTT        0       40690558        T       A
# 13      chr13:40690558:A:ATTTTTTT       0       40690558        T       A
# 13      chr13:40690558:AT:A     0       40690558        A       T

plink --bfile prs_out_wgs_survivor/POI_META_harmonized --freq -out prs_out_wgs_survivor/POI_META_harmonized_freq 
## This is a multiallelic variant, I am keeping only this variant: chr13:40690558:AT:A
 # CHR                         SNP   A1   A2          MAF  NCHROBS
 #  13         chr13:40690558:A:AT    T    A      0.04103    15938
 #  13      chr13:40690558:A:ATTTT    T    A      0.01449    15938
 #  13     chr13:40690558:A:ATTTTT    T    A      0.07059    15938
 #  13    chr13:40690558:A:ATTTTTT    T    A       0.1996    15938
 #  13   chr13:40690558:A:ATTTTTTT    T    A       0.0155    15938
 #  13         chr13:40690558:AT:A    A    T       0.1579    15938

echo "chr13:40690558:AT:A" > prs_out_wgs_survivor/${study}_snp_list.txt
plink --bfile prs_out_wgs_survivor/POI_META_harmonized_tmp  --extract prs_out_wgs_survivor/${study}_snp_list.txt --make-bed --out prs_out_wgs_survivor/${study}_harmonized ##**

## duplicate variants
plink --bfile prs_out_wgs_survivor/${study}_direct_match --bmerge prs_out_wgs_survivor/${study}_harmonized --make-bed --out prs_out_wgs_survivor/$study ##**

## Check if any duplicates
awk '{print $1 ":" $4}' prs_out_wgs_survivor/$study.bim | sort | uniq -c | sort -nr | head
# (base) [aneupane@splprhpc12 prs]$ awk '{print $1 ":" $4}' prs_out_wgs_survivor/$study.bim | sort | uniq -c | sort -nr | head
#       2 9:134097695
#       2 7:55967632
#       2 7:5405116
#       2 6:158153983
#       2 14:103553662
#       1 9:878563
#       1 9:33012384
#       1 9:23823401
#       1 9:123797714
#       1 8:94619248
# (base) [aneupane@splprhpc12 prs]$ grep 134097695 prs_out_wgs_survivor/$study.bim
# 9       chr9:134097695:C:CA     0       134097695       CA      C
# 9       chr9:134097695:CA:C     0       134097695       C       CA
# (base) [aneupane@splprhpc12 prs]$ grep 55967632 prs_out_wgs_survivor/$study.bim
# 7       chr7:55967632:A:ATTAT   0       55967632        ATTAT   A
# 7       chr7:55967632:ATTAT:A   0       55967632        A       ATTAT
# (base) [aneupane@splprhpc12 prs]$ grep 5405116 prs_out_wgs_survivor/$study.bim
# 7       chr7:5405116:C:CA       0       5405116 CA      C
# 7       chr7:5405116:CA:C       0       5405116 C       CA
# (base) [aneupane@splprhpc12 prs]$ grep 158153983 prs_out_wgs_survivor/$study.bim
# 6       chr6:158153983:A:AAAT   0       158153983       AAAT    A
# 6       chr6:158153983:AAAT:A   0       158153983       A       AAAT
# (base) [aneupane@splprhpc12 prs]$ grep 103553662 prs_out_wgs_survivor/$study.bim
# 14      chr14:103553662:G:GA    0       103553662       GA      G
# 14      chr14:103553662:GA:G    0       103553662       G       GA


plink --bfile prs_out_wgs_survivor/$study --freq --out prs_out_wgs_survivor/$study_freq

# remove these:
# chr9:134097695:CA:C    C   CA     0.005767    15954
# chr7:55967632:ATTAT:A    A ATTAT        0.127    15954
# chr7:5405116:C:CA   CA    C      0.03642    15954
# chr6:158153983:A:AAAT AAAT    A    0.0003134    15954
# chr14:103553662:GA:G    G   GA      0.00539    15954

cat <<EOT > prs_out_wgs_survivor/remove_variants.txt
chr9:134097695:CA:C
chr7:55967632:ATTAT:A
chr7:5405116:C:CA
chr6:158153983:A:AAAT
chr14:103553662:GA:G
EOT

plink --bfile prs_out_wgs_survivor/${study} --exclude  prs_out_wgs_survivor/remove_variants.txt --make-bed --out prs_out_wgs_survivor/${study}2

# Update variant names
awk '{print $2, "chr"$1":"$4}' prs_out_wgs_survivor/${study}2.bim > prs_out_wgs_survivor/${study}_update_variantnames
plink --bfile prs_out_wgs_survivor/${study}2 --update-name prs_out_wgs_survivor/${study}_update_variantnames --make-bed --out prs_out_wgs_survivor/${study}_varname_updated



# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out_wgs_survivor/POI_meta_GRCh38.dat_${study} > prs_out_wgs_survivor/${study}.prsweight

# Calculate PRS
plink --bfile prs_out_wgs_survivor/${study}_varname_updated --score prs_out_wgs_survivor/${study}.prsweight --out prs_out_wgs_survivor/${study}_prs
# --score: 279 valid predictors loaded.





#####################################
## create plink files for CCSS org ##
#####################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/plink_data/

module load plink/1.90b
for CHR in {1..22} X Y; do
echo "Doing Chr ${CHR}"
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/CCSS_org_GRCh37_chr${CHR} --extract ../wantedSNP_GRCh37.txt --keep-allele-order --make-bed --out Ruth_etal_PRS_chr${CHR}_ccss_org
done


ls Ruth_etal_PRS_chr*_ccss_org.bed | sed 's/\..*//' | sort -u > merge_list.txt
plink --merge-list merge_list.txt -keep-allele-order --make-bed --out Ruth_etal_PRS_merged_ccss_org


## After running R code PRS_scores.R; extract variants
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/
ln -s  ../POI_meta_GRCh38.dat .
ln -s  ../POI_meta_GRCh37.dat .

ln -s ../plink_data/Ruth_etal_PRS_merged_ccss_org.* .

## After running R code PRS_scores.R; extract variants

## PRS; studies: PGS000356 (172 vars out of 179) PGS000454 (27 vars out of 27) PGS003416 (460 vars out of 462) ST6 (110 vars out of 110)
mkdir -p prs_out_ccss_org
study="POI_META"
# Subset PRS data for each study
awk -v study=$study '$6==study' POI_meta_GRCh37.dat > prs_out_ccss_org/POI_meta_GRCh37.dat_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out_ccss_org/POI_meta_GRCh37.dat_${study} | wc -l
## Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out_ccss_org/POI_meta_GRCh37.dat_${study} Ruth_etal_PRS_merged_ccss_org.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_direct_match

# (base) [aneupane@splprhpc12 prs]$ wc -l prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_direct_match
# 283 prs_out_ccss_org/POI_meta_GRCh37.dat_POI_META_direct_match

# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out_ccss_org/POI_meta_GRCh37.dat_${study} Ruth_etal_PRS_merged_ccss_org.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v '*' > prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_direct_match prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match \
> prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match_final

wc -l prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match_final
# 6 prs_out_ccss_org/POI_meta_GRCh37.dat_POI_META_no_direct_match_final

## Harmonize alleles
module load R/4.2.2-rhel8
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match_final

## Update the alleles; It checks if the last field ($NF) in each line is equal to 1; if true, then { print $3, $6, $7, $8, $9 } so we can update the plink file for variants with no direct match
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match_final_alleles_harmonized > prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt


# Extract study-specific variants
awk '{print $2}' prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_direct_match > prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_direct_match_to_extract.txt
module load plink/1.90b


module load plink/1.90b

## If there are no indirect match, just make this commented first line: --out prs_out_ccss_org/${study} and skip lines with **
# plink --bfile Ruth_etal_PRS_merged_ccss_org --extract prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out_ccss_org/${study}_direct_match
plink --bfile Ruth_etal_PRS_merged_ccss_org --extract prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out_ccss_org/${study}
# plink --bfile Ruth_etal_PRS_merged_ccss_org --extract prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out_ccss_org/POI_meta_GRCh37.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out_ccss_org/${study}_harmonized ##**
# This command above extracts and updates alleles from 
# (base) [aneupane@splprhpc08 prs_out_ccss_org]$ grep chr10:46037697 sjlife_extracted_chrALL.bim
# 10      chr10:46037697:A:G      0       46037697        G       A
# to ---------->
# (base) [aneupane@splprhpc08 prs_out_ccss_org]$ cat Pleiotropy_One_directional_Significant_harmonized.bim
# 10      chr10:46037697:A:G      0       46037697        C       T

# plink --bfile prs_out_ccss_org/${study}_direct_match --bmerge prs_out_ccss_org/${study}_harmonized --make-bed --out prs_out_ccss_org/$study ##**
# Update variant names
awk '{print $2, "chr"$1":"$4}' prs_out_ccss_org/${study}.bim > prs_out_ccss_org/${study}_update_variantnames
# wc -l prs_out_ccss_org/${study}_update_variantnames
# 231 prs_out_ccss_org/POI_META_update_variantnames

plink --bfile prs_out_ccss_org/$study --update-name prs_out_ccss_org/${study}_update_variantnames --make-bed --out prs_out_ccss_org/${study}_varname_updated
# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out_ccss_org/POI_meta_GRCh37.dat_${study} > prs_out_ccss_org/${study}.prsweight

# Calculate PRS
plink --bfile prs_out_ccss_org/${study}_varname_updated --score prs_out_ccss_org/${study}.prsweight --out prs_out_ccss_org/${study}_prs


