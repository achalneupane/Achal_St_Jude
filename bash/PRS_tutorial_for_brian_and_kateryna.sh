## Step 1.
# extract GWAS variants from Michighanweb server (with the highest AUC) or PRScatalog 

# step 2,  
# GRCh37 ---> GRCh38

# Step 3. 
# Extract variants from your data

## Fix effect size; has to be positive effect size consistently
head(all.cancer)
#   CHROM POS_GRCh38 REF Effect_allele Effect_size                                 TYPE     Cancer Significant_YN
# 1 chr10  121580917   G            GC  0.24302461 Pleiotropy_Replication_prior_studies Pleiotropy              Y
# 2 chr10  121566241   C             T -0.06315282 Pleiotropy_Replication_prior_studies Pleiotropy              Y
# 3 chr16   52477173   A             G -0.22765371 Pleiotropy_Replication_prior_studies Pleiotropy              Y
# 4 chr16   52528555   T             G -0.21828040 Pleiotropy_Replication_prior_studies Pleiotropy              Y
# 5 chr16   52565276   T             C -0.23787655 Pleiotropy_Replication_prior_studies Pleiotropy              Y

# <R>
all.cancers[all.cancers$Effect_size < 0, c("REF", "Effect_allele")] <- all.cancers[all.cancers$Effect_size < 0, c("Effect_allele", "REF")]
all.cancers$Effect_size <- abs(all.cancers$Effect_size)

# step 4. 
## Create bed files from all.cancer
# (base) [aneupane@splprhpc08 plink_data]$ head PRS_vars.bed
# chr11   69516649        69516650
# chr11   69564392        69564393
# chr11   69516873        69516874
# chr5    1282203 1282204
# chr5    1279674 1279675
# chr5    1297372 1297373

# <sh>
hpcf_interactive -n 12 -R "rusage[mem=8000]" -q standard

#!/usr/bin/bash
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/plink_data

# load modules
module load bcftools/1.9
module load plink/1.90b


THREADS=4
for CHR in {1..22}; do
grep -w chr${CHR} PRS_vars.bed | sort -V > PRS_vars_chr${CHR}.bed
sed -i "s/\r//g" PRS_vars_chr${CHR}.bed
bcftools view -Oz MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz --threads ${THREADS} -R PRS_vars_chr${CHR}.bed > PRS_chr${CHR}.vcf.gz
plink --vcf PRS_chr${CHR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out PRS_chr${CHR}
done


## step 6
## Harmonize alleles and run PRS calculation
ALL_Cancers_PRS_data.txt ## you only need CHR, POS, Effect allele and Effect size.
(base) [aneupane@splprhpc08 prs]$ head ALL_Cancers_PRS_data.txt
# CHROM   POS_GRCh38      REF     Effect_allele   Effect_size     TYPE    Cancer  Significant_YN
# 11      69516650        T       C       0.099930839     Mavaddat_2015_ER_NEG_Breast     Breast  Y
# 11      69564393        C       A       0.124515875     Mavaddat_2015_ER_NEG_Breast     Breast  Y
# 11      69516874        C       G       0.009356095     Mavaddat_2015_ER_NEG_Breast     Breast  Y
# 5       1282204 C       A       0.030141157     Mavaddat_2015_ER_NEG_Breast     Breast  Y
# 5       1279675 C       T       0.113328685     Mavaddat_2015_ER_NEG_Breast     Breast  Y
# 5       1297373 T       C       0.109368537     Mavaddat_2015_ER_NEG_Breast     Breast  Y

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs
mkdir -p prs_out
study="Mavaddat_2019_ER_NEG_Breast"
# Subset PRS data for each study
awk -v study=$study '$6==study' ALL_Cancers_PRS_data.txt > prs_out/ALL_Cancers_PRS_data.txt_${study}
# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/ALL_Cancers_PRS_data.txt_$study | wc -l
## Look for directly matching variants in the WGS data
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/ALL_Cancers_PRS_data.txt_$study plink_data/sjlife_all_PRS_all_final_v3.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match

## Code description from above:
# (NR==FNR) checks whether the total number of records processed (NR) is equal to the number of records processed in the current input file (FNR).
# NR==FNR{a[$1":"$2]=$3" "$4;next}:
# While reading the first file (NR==FNR), it creates a memory storage (a) where it combines the first and second columns of the file (e.g., chr1:12345) as the key and stores the third and fourth columns (e.g., A T) as the value.
# The next command tells AWK to move to the next line and not process anything else in the code for this file.
# ($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}:
# For the second file (when NR is not equal to FNR), it checks if the combination of the first and fourth columns (e.g., chr1:67890) exists as a key in the a memory storage created earlier.
# If it finds a match, it prints the first, second, fourth, fifth, and sixth columns from the current line, along with the corresponding value from memory (a[$1":"$4]) which is the third and fourth columns from the first file (e.g., A T).
# Essentially, it combines information from the two files based on matching keys from the first and fourth columns.

## This part takes the output of the first awk command and performs additional filtering:
# | awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)' > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match
# ($4==$6 || $4==$7) && ($5==$6 || $5==$7): It checks if either the fourth column ($4) matches either the sixth column ($6) or the seventh column ($7), and if the fifth column ($5) also matches either the sixth column ($6) or the seventh column ($7).
# If this condition is met, it passes the line to the output.
# Finally, the output is redirected to a file named prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match.


# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/ALL_Cancers_PRS_data.txt_$study plink_data/sjlife_all_PRS_all_final_v3.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match \
> prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

wc -l prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

## Harmonize alleles
module load R
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final

## Update the alleles; It checks if the last field ($NF) in each line is equal to 1; if true, then { print $3, $6, $7, $8, $9 } so we can update the plink file for variants with no direct match
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_final_alleles_harmonized > prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# (base) [aneupane@splprhpc08 prs_out]$ cat ALL_Cancers_PRS_data.txt_Pleiotropy_One_directional_Significant_no_direct_match_final_alleles_harmonized
# chr pos variant khera_a1 khera_a2 wgs_a1 wgs_a2 wgs_a1_new wgs_a2_new match
# 10 46037697 chr10:46037697:A:G T C G A C T 1

# (base) [aneupane@splprhpc08 prs_out]$ cat ALL_Cancers_PRS_data.txt_Pleiotropy_One_directional_no_direct_match_alleles_harmonized_update_alleles.txt
# chr10:46037697:A:G G A C T


# Extract study-specific variants
awk '{print $2}' prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match > prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match_to_extract.txt
module load plink/1.90b
plink --bfile plink_data/sjlife_all_PRS_all_final_v3 --extract prs_out/ALL_Cancers_PRS_data.txt_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match
plink --bfile plink_data/sjlife_all_PRS_all_final_v3 --extract prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/ALL_Cancers_PRS_data.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized
# This command above extracts and updates alleles from 
# (base) [aneupane@splprhpc08 prs_out]$ grep chr10:46037697 ../plink_data/sjlife_all_PRS_all_final_v3.bim
# 10      chr10:46037697:A:G      0       46037697        G       A
# to ---------->
# (base) [aneupane@splprhpc08 prs_out]$ cat Pleiotropy_One_directional_Significant_harmonized.bim
# 10      chr10:46037697:A:G      0       46037697        C       T

plink --bfile prs_out/${study}_direct_match --bmerge prs_out/${study}_harmonized --make-bed --out prs_out/$study
# Update variant names
awk '{print $2, $1":"$4}' prs_out/${study}.bim > prs_out/${study}_update_variantnames
plink --bfile prs_out/$study --update-name prs_out/${study}_update_variantnames --make-bed --out prs_out/${study}_varname_updated
# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out/ALL_Cancers_PRS_data.txt_$study > prs_out/${study}.prsweight
# (base) [aneupane@splprhpc08 prs_out]$ head Pleiotropy_One_directional.prsweight
# 8:127401060 G 0.214304603
# 17:37743574 A 0.196014884
# 19:50862184 C 0.337899789
# 6:32623790 T 0.392717535
# 11:69213892 G 0.223143551


# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs

# Done!