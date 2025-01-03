# Dear Achal,

# I am writing to ask you for your help. Do you think you might be able to help us put the following together over the next 2 weeks to meet the mid-Feb ISLCCC abstract deadline?

# Basically, we would like to get P/LP variant carrier variables in SJLIFE (WES/WGS), CCSS Expansion (WES/WGS), and CCSS Original (WES) considering the following 5 gene lists below, using any of the 3 rare variant masks that were described in our POI concept proposal (SnpEff; LOFTEE; ClinVar - I think if we have to pick one, we should pick the strictest definition, which is probably ClinVar). Since this is a quick preliminary pass, maybe we could focus on getting at least one cancer susceptibility P/LP variant carrier status variable ([1] or [2] below) and at least one BCC-related variable ([3]-[5] below). Do you think this would be feasible? Please let me know if you have questions or if you would like to meet to discuss.

# Thank you so much,
# Cindy

# Gene lists 

# (1) 60 cancer susceptibility genes - should be Zhaoming's SJCPG60 list. See Kim_ST1.txt, use "CSG_60" (genes marked with "x"). From Kim et al evaluating cancer susceptibility gene P/LP variants in CCSS, link: https://doi.org/10.1093/jncics/pkab007.

# (2) Expanded list of 172 cancer susceptibility genes. Kim_ST1.txt, use "CSG_172" (genes marked with "x"; 172 genes evaluated in CCSS)

# (3) Literature-based list of BCC genes. See "NCI table 3" for a list of ~20 genes (Table 3 in the cancer.gov). These genes are based on what has been reported in the literature re: BCC-associated syndromes.
# Based on Kilgour et al, Choquet et al, and NCI cancer.gov (all 3 are very consistent):
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8345475/#B6-cancers-13-03870
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7259534/
# https://www.cancer.gov/types/skin/hp/skin-genetics-pdq

# (4) BCC gene panel genes. See "BCC_blueprint_panel.txt" for a list of genes tested by Blueprint Genetics panel for BCC:
# https://blueprintgenetics.com/tests/panels/dermatology/hereditary-melanoma-and-skin-cancer-panel/

# (5) ClinVar BCC-related genes. See "clinvar_result_BCC.txt". These are all ClinVar listings for "basal AND cell AND carcinoma" annotated as P or LP and with multiple submitters/no conflicts. This list hasn't been reviewed yet in terms of conditions, but I think we can use the variant/gene list. For the final analysis, I can ask the clinicians to help review the conditions to make sure we don't have anything unrelated.


# Further details in Email from 1/25/2024


## Extract variants for PRS
# load modules
module load bcftools
module load plink/1.90b


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/prs/sjlife
###################################
## create plink files for sjlife ##
###################################
THREADS=4
for CHR in {1..22}; do
grep -w chr${CHR} ../all_bed_BCC.bed | sort -V > PRS_vars_chr${CHR}.bed
sed -i "s/\r//g" PRS_vars_chr${CHR}.bed
bcftools view -Oz /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2//MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz --threads ${THREADS} -R PRS_vars_chr${CHR}.bed > PRS_chr${CHR}.vcf.gz
plink --vcf PRS_chr${CHR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out PRS_chr${CHR}
done


## After running R code PRS_scores.R; extract variants
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/prs/sjlife
ln -s ../all_BCC_effect.dat .
ln -s ../extract_sjlife_vars .


## using preQC data
for CHR in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz --extract extract_sjlife_vars --keep-allele-order --make-bed --out sjlife_extracted_chr${CHR}
done

ls sjlife_extracted_chr*.bim | sed 's/.bim//'| sort -V > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out sjlife_extracted_chrALL


# Append chr in the first column
awk '{print "chr"$1, $2, $3, $4, $5, $6}' sjlife_extracted_chrALL.bim > tmp_bim
mv tmp_bim sjlife_extracted_chrALL.bim

## Note: SJLIFE and CCSS_exp have all variants; CCSS_org is missing 9 variants in total
## PRS; studies: PGS000356 (179 vars) PGS000454 (27 vars) PGS003416 (462 vars) ST6 (110 vars)
mkdir -p prs_out
study="PGS003416"
# Subset PRS data for each study
awk -v study=$study '$6==study' all_BCC_effect.dat > prs_out/all_BCC_effect.dat_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/all_BCC_effect.dat_${study} | wc -l
## Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/all_BCC_effect.dat_${study} sjlife_extracted_chrALL.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/all_BCC_effect.dat_${study}_direct_match


# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/all_BCC_effect.dat_${study} sjlife_extracted_chrALL.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/all_BCC_effect.dat_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/all_BCC_effect.dat_${study}_direct_match prs_out/all_BCC_effect.dat_${study}_no_direct_match \
> prs_out/all_BCC_effect.dat_${study}_no_direct_match_final

wc -l prs_out/all_BCC_effect.dat_${study}_no_direct_match_final

## Harmonize alleles
module load R
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R prs_out/all_BCC_effect.dat_${study}_no_direct_match_final

## Update the alleles; It checks if the last field ($NF) in each line is equal to 1; if true, then { print $3, $6, $7, $8, $9 } so we can update the plink file for variants with no direct match
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/all_BCC_effect.dat_${study}_no_direct_match_final_alleles_harmonized > prs_out/all_BCC_effect.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt


# Extract study-specific variants
awk '{print $2}' prs_out/all_BCC_effect.dat_${study}_direct_match > prs_out/all_BCC_effect.dat_${study}_direct_match_to_extract.txt
module load plink/1.90b

## If there are no indirect match, just make this commented first line: --out prs_out/${study} and skip lines with **
# plink --bfile sjlife_extracted_chrALL --extract prs_out/all_BCC_effect.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
plink --bfile sjlife_extracted_chrALL --extract prs_out/all_BCC_effect.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}
plink --bfile sjlife_extracted_chrALL --extract prs_out/all_BCC_effect.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/all_BCC_effect.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized ##**
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
awk '{print $1":"$2, $4, $5}' prs_out/all_BCC_effect.dat_${study} > prs_out/${study}.prsweight

# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs





#####################################
## create plink files for CCSS exp ##
#####################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome
for i in {1..22}; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
unset VCF; \
export THREADS=4; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/"; \
bsub \
        -P "${CHR}_plk" \
        -J "${CHR}_plk" \
        -o "${WORKDIR}/logs/${VCF%.vcf*}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "./plinkQC.sh"; \
done;

#!/usr/bin/bash
module load plink/1.90b
plink --vcf "${WORKDIR}/CCSS_exp_biallelic_${CHR}_ID_edited.vcf.gz" --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out "${WORKDIR}/CCSS_exp_biallelic_${CHR}_ID_edited"




## After running R code PRS_scores.R; extract variants
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_exp
ln -s ../all_BCC_effect.dat .
ln -s ../extract_CCSS_exp_vars .



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
study="ST6"
# Subset PRS data for each study
awk -v study=$study '$6==study' all_BCC_effect.dat > prs_out/all_BCC_effect.dat_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/all_BCC_effect.dat_${study} | wc -l
## Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/all_BCC_effect.dat_${study} ccss_exp_extracted_chrALL.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/all_BCC_effect.dat_${study}_direct_match


# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/all_BCC_effect.dat_${study} ccss_exp_extracted_chrALL.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/all_BCC_effect.dat_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/all_BCC_effect.dat_${study}_direct_match prs_out/all_BCC_effect.dat_${study}_no_direct_match \
> prs_out/all_BCC_effect.dat_${study}_no_direct_match_final

wc -l prs_out/all_BCC_effect.dat_${study}_no_direct_match_final

## Harmonize alleles
module load R
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R prs_out/all_BCC_effect.dat_${study}_no_direct_match_final

## Update the alleles; It checks if the last field ($NF) in each line is equal to 1; if true, then { print $3, $6, $7, $8, $9 } so we can update the plink file for variants with no direct match
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/all_BCC_effect.dat_${study}_no_direct_match_final_alleles_harmonized > prs_out/all_BCC_effect.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt


# Extract study-specific variants
awk '{print $2}' prs_out/all_BCC_effect.dat_${study}_direct_match > prs_out/all_BCC_effect.dat_${study}_direct_match_to_extract.txt
module load plink/1.90b

## If there are no indirect match, just make this commented first line: --out prs_out/${study} and skip lines with **
# plink --bfile ccss_exp_extracted_chrALL --extract prs_out/all_BCC_effect.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
plink --bfile ccss_exp_extracted_chrALL --extract prs_out/all_BCC_effect.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}
plink --bfile ccss_exp_extracted_chrALL --extract prs_out/all_BCC_effect.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/all_BCC_effect.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized ##**
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
awk '{print $1":"$2, $4, $5}' prs_out/all_BCC_effect.dat_${study} > prs_out/${study}.prsweight

# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs




#####################################
## create plink files for CCSS org ##
#####################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/
for i in {1..22}; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
export THREADS=4; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/"; \
export VCFDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/"; \
bsub \
        -P "${CHR}_plk" \
        -J "${CHR}_plk" \
        -o "${WORKDIR}/logs/${CHR}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "./plinkQC.sh"; \
done;

# <plinkQC.sh>
#!/usr/bin/bash
module load plink/1.90b
plink --vcf "${VCFDIR}/${CHR}.dose.vcf.gz" --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out "${WORKDIR}/CCSS_org_GRCh37_${CHR}"



# #!/usr/bin/bash
# module load bcftools/1.9
# module load tabix
# cd ${WORKDIR}
# /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/liftover/liftOver "${VCFDIR}/${CHR}.dose.vcf.gz" /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/liftover/hg19ToHg38.over.chain ${CHR}_GRCh38.dose.vcf.gz ${CHR}_unmapped.bed
# tabix -p vcf "${WORKDIR}/${CHR}_GRCh38.dose.vcf.gz"


## Rename all BIM files with CHR:POS:REF:ALT in shell script in ccss_org
for i in {1..22}; do \
CHR="chr${i}";
awk '{print "chr"$1, "chr"$2":"$6":"$5, $3, $4, $5, $6}' CCSS_org_GRCh37_${CHR}.bim > tmp_bim
mv tmp_bim CCSS_org_GRCh37_${CHR}.bim
done


## After running R code PRS_scores.R; extract variants
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_org
ln -s ../all_BCC_effect_ccss_org.dat .
ln -s ../extract_ccss_org_vars .



## using preQC data
for CHR in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/CCSS_org_GRCh37_chr${CHR} --extract extract_ccss_org_vars --keep-allele-order --make-bed --out ccss_org_extracted_chr${CHR}
done


ls ccss_org_extracted_chr*.bim | sed 's/.bim//'| sort -V > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out ccss_org_extracted_chrALL

# Append chr in the first column
awk '{print "chr"$1, $2, $3, $4, $5, $6}' ccss_org_extracted_chrALL.bim > tmp_bim
mv tmp_bim ccss_org_extracted_chrALL.bim

## PRS; studies: PGS000356 (172 vars out of 179) PGS000454 (27 vars out of 27) PGS003416 (460 vars out of 462) ST6 (110 vars out of 110)
mkdir -p prs_out
study="ST6"
# Subset PRS data for each study
awk -v study=$study '$6==study' all_BCC_effect_ccss_org.dat > prs_out/all_BCC_effect_ccss_org.dat_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/all_BCC_effect_ccss_org.dat_${study} | wc -l
## Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/all_BCC_effect_ccss_org.dat_${study} ccss_org_extracted_chrALL.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/all_BCC_effect_ccss_org.dat_${study}_direct_match


# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/all_BCC_effect_ccss_org.dat_${study} ccss_org_extracted_chrALL.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/all_BCC_effect_ccss_org.dat_${study}_direct_match prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match \
> prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match_final

wc -l prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match_final

## Harmonize alleles
module load R
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match_final

## Update the alleles; It checks if the last field ($NF) in each line is equal to 1; if true, then { print $3, $6, $7, $8, $9 } so we can update the plink file for variants with no direct match
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match_final_alleles_harmonized > prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt


# Extract study-specific variants
awk '{print $2}' prs_out/all_BCC_effect_ccss_org.dat_${study}_direct_match > prs_out/all_BCC_effect_ccss_org.dat_${study}_direct_match_to_extract.txt
module load plink/1.90b

## If there are no indirect match, just make this commented first line: --out prs_out/${study} and skip lines with **
# plink --bfile ccss_org_extracted_chrALL --extract prs_out/all_BCC_effect_ccss_org.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
plink --bfile ccss_org_extracted_chrALL --extract prs_out/all_BCC_effect_ccss_org.dat_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}
plink --bfile ccss_org_extracted_chrALL --extract prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/all_BCC_effect_ccss_org.dat_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized ##**
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
awk '{print $1":"$2, $4, $5}' prs_out/all_BCC_effect_ccss_org.dat_${study} > prs_out/${study}.prsweight

# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs


