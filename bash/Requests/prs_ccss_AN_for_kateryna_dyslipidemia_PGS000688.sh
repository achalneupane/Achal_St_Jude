##############
## CCSS_exp ##
##############

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_exp
## using preQC data
for CHR in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited --extract extract_CCSS_exp_vars --keep-allele-order --make-bed --out ccss_exp_extracted_chr${CHR}
done


ls ccss_exp_extracted_chr*.bim | sed 's/.bim//'| sort -V > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out ccss_exp_extracted_chrALL

# Append chr in the first column
awk '{print "chr"$1, $2, $3, $4, $5, $6}' ccss_exp_extracted_chrALL.bim > tmp_bim
mv tmp_bim ccss_exp_extracted_chrALL.bim

## PRS; studies: PGS000356 (179 vars) PGS000454 (27 vars) PGS003416 (462 vars) ST6 (110 vars)
mkdir -p prs_out
study="PGS000688"
# Subset PRS data for each study
awk -v study=$study '$6==study' GRCh38_PGS000688.dat > prs_out/GRCh38_PGS000688.dat_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/GRCh38_PGS000688.dat_${study} | wc -l
## Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/GRCh38_PGS000688.dat_${study} ccss_exp_extracted_chrALL.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/GRCh38_PGS000688.dat_${study}_direct_match


# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/GRCh38_PGS000688.dat_${study} ccss_exp_extracted_chrALL.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/GRCh38_PGS000688.dat_${study}_direct_match prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match \
> prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match_final

wc -l prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match_final

## Harmonize alleles
module load R/4.0.2-a-rhel8
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match_final

## Update the alleles; It checks if the last field ($NF) in each line is equal to 1; if true, then { print $3, $6, $7, $8, $9 } so we can update the plink file for variants with no direct match
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match_final_alleles_harmonized > prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt


# Extract study-specific variants
awk '{print $2}' prs_out/GRCh38_PGS000688.dat_${study}_direct_match > prs_out/GRCh38_PGS000688.dat_${study}_direct_match_to_extract.txt
module load plink/1.90b

## If there are no indirect match, just make this commented first line: --out prs_out/${study} and skip lines with **
plink --bfile ccss_exp_extracted_chrALL --extract prs_out/GRCh38_PGS000688.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
# plink --bfile ccss_exp_extracted_chrALL --extract prs_out/GRCh38_PGS000688.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}
plink --bfile ccss_exp_extracted_chrALL --extract prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/GRCh38_PGS000688.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized ##**


plink --bfile prs_out/${study}_direct_match --bmerge prs_out/${study}_harmonized --make-bed --out prs_out/$study ##**
# Update variant names
awk '{print $2, "chr"$1":"$4}' prs_out/${study}.bim > prs_out/${study}_update_variantnames
plink --bfile prs_out/$study --update-name prs_out/${study}_update_variantnames --make-bed --out prs_out/${study}_varname_updated
# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out/GRCh38_PGS000688.dat_${study} > prs_out/${study}.prsweight

# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs

##############
## CCSS_org ##
##############
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_org
## using preQC data
for CHR in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/CCSS_org_GRCh37_chr${CHR} --extract extract_CCSS_org_vars --keep-allele-order --make-bed --out ccss_org_extracted_chr${CHR}
done


ls ccss_org_extracted_chr*.bim | sed 's/.bim//'| sort -V > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out ccss_org_extracted_chrALL

# Append chr in the first column
awk '{print "chr"$1, $2, $3, $4, $5, $6}' ccss_org_extracted_chrALL.bim > tmp_bim
mv tmp_bim ccss_org_extracted_chrALL.bim

## PRS; studies: PGS000356 (172 vars out of 179) PGS000454 (27 vars out of 27) PGS003416 (460 vars out of 462) ST6 (110 vars out of 110)
mkdir -p prs_out
study="PGS000688"
# Subset PRS data for each study
awk -v study=$study '$6==study' GRCh37_PGS000688.dat > prs_out/GRCh37_PGS000688.dat_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/GRCh37_PGS000688.dat_${study} | wc -l
## Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/GRCh37_PGS000688.dat_${study} ccss_org_extracted_chrALL.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/GRCh37_PGS000688.dat_${study}_direct_match


# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/GRCh37_PGS000688.dat_${study} ccss_org_extracted_chrALL.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/GRCh37_PGS000688.dat_${study}_direct_match prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match \
> prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match_final

wc -l prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match_final

## Harmonize alleles
module load R
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match_final

## Update the alleles; It checks if the last field ($NF) in each line is equal to 1; if true, then { print $3, $6, $7, $8, $9 } so we can update the plink file for variants with no direct match
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match_final_alleles_harmonized > prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt


# Extract study-specific variants
awk '{print $2}' prs_out/GRCh37_PGS000688.dat_${study}_direct_match > prs_out/GRCh37_PGS000688.dat_${study}_direct_match_to_extract.txt
module load plink/1.90b

## If there are no indirect match, just make this commented first line: --out prs_out/${study} and skip lines with **
# plink --bfile ccss_org_extracted_chrALL --extract prs_out/GRCh37_PGS000688.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
plink --bfile ccss_org_extracted_chrALL --extract prs_out/GRCh37_PGS000688.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}
plink --bfile ccss_org_extracted_chrALL --extract prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/GRCh37_PGS000688.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized ##**
# This command above extracts and updates alleles from 
# (base) [aneupane@splprhpc08 prs_out]$ grep chr10:46037697 sjlife_extracted_chrALL.bim
# 10      chr10:46037697:A:G      0       46037697        G       A
# to ---------->
# (base) [aneupane@splprhpc08 prs_out]$ cat Pleiotropy_One_directional_Significant_harmonized.bim
# 10      chr10:46037697:A:G      0       46037697        C       T

plink --bfile prs_out/${study}_direct_match --bmerge prs_out/${study}_harmonized --make-bed --out prs_out/$study ##**
# Update variant names
awk '{print $2, "chr"$1":"$4}' prs_out/${study}.bim > prs_out/${study}_update_variantnames
plink --bfile prs_out/$study --update-name prs_out/${study}_update_variantnames --make-bed --out prs_out/${study}_varname_updated
# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out/GRCh37_PGS000688.dat_${study} > prs_out/${study}.prsweight

# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs