####################
## CCSS_Expansion ##
####################
## extract variants from VCF and see which ones are present
module load bcftools/1.9
module load plink/1.90b


# split chromosomes
for i in {1..22}; do \
chr=$(echo "chr"${i}); \
export CHR=${chr}; \
echo $CHR
export MEM=6; \
export THREADS=4; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/"; \
export OUTDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/"; \
bsub \
        -P "chr${CHR}_extract" \
        -J "chr${CHR}_extract" \
        -o "${OUTDIR}/logs/${CHR}_s00VCFextract.%J" \
        -n ${THREADS} \
        -R "rusage[mem=6192]" \
        "./extract_chr.sh"; \
done;

cat <<\EoF > extract_chr.sh
 #!/usr/bin/bash
module load bcftools/1.9
module load tabix
bcftools view ${WORKDIR}/CCSS.vcf.gz --regions ${CHR} | bgzip -c > ${OUTDIR}/splitted_${CHR}_CCSS.vcf.gz
tabix -p vcf ${OUTDIR}/splitted_${CHR}_CCSS.vcf.gz
EoF




# Make VCF biallelic and also edit SNP IDs and header
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/
cat <<\EoF > edit_vcf_header.sh
#!/usr/bin/bash
module load bcftools/1.10.2

cd "${OUT_DIR}"
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa splitted_${CHR}_CCSS.vcf.gz -Oz -o CCSS_exp_biallelic_${CHR}_tmp1.vcf.gz
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' CCSS_exp_biallelic_${CHR}_tmp1.vcf.gz -Oz -o CCSS_exp_biallelic_${CHR}_ID_edited_tmp2.vcf.gz
bcftools reheader -s COMPBIO_ID_TO_CCSS_ID.txt CCSS_exp_biallelic_${CHR}_ID_edited_tmp2.vcf.gz > CCSS_exp_biallelic_${CHR}_ID_edited.vcf.gz
bcftools index -f -t --threads 4 CCSS_exp_biallelic_${CHR}_ID_edited.vcf.gz
EoF



for i in {1..22}; do \
chr=$(echo "chr"${i}); \
export CHR=${chr}; \
echo $CHR; \
export MEM=6; \
export THREADS=4; \
export OUTDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/"; \
bsub \
        -P "chr${CHR}_extract" \
        -J "chr${CHR}_extract" \
        -o "${OUTDIR}/logs/${CHR}_header.%J" \
        -n ${THREADS} \
        -R "rusage[mem=6192]" \
        "./edit_vcf_header.sh"; \
done;



# extract variants from preQC VCF
# first, run extract_variants.py to get variant and PRS score. Then:
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs

module load bcftools/1.9
module load plink/1.90b
for line in $(cat /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/all_cancer_extract_var.txt); do
VAR="$(echo ${line}| tr -d " \t\n\r" )"
CHR="$(echo $VAR |awk -F':' '{print $1}')"
echo "Doing ${VAR}"
bcftools view /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_${CHR}_ID_edited.vcf.gz ${VAR} > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR}.vcf.gz
plink --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --make-bed --out /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR} 2>&1 | tee -a extract_plink_all.log
done

cd plink_data
# find which ones are missing
for line in $(cat ../all_cancer_extract_var.txt); do
VAR="$(echo ${line}| tr -d " \t\n\r" )"
CHR="$(echo $VAR |awk -F':' '{print $1}')"
if [ ! -f "PRS_${VAR}.bim" ] ; then
echo ${VAR} >> failed_plink_files.txt
fi
done

ls *.bim| sort -V | sed 's/\.bim//g'|sed -n '1d;p' > merge_list.list
plink --bfile PRS_chr1:2125052 --merge-list merge_list.list --keep-allele-order --out merged.dat

# ------------------------------
# Now run process_prs.sh >>
#!/bin/bash
# ------------------------------
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs

mkdir prs_out

study=$1

# ALL_Vijayakrishnan
# Basal_cell_carcinoma_PRSWeb
# Squamous_cell_carcinoma_PRSWeb
# Mavaddat_2019_ER_NEG_Breast
# Mavaddat_2019_ER_OVERALL_Breast
# Mavaddat_2019_ER_POS_Breast
# Meningioma
# Pleiotropy_PRSWEB
# Sarcoma_Machiela
# THYROID_PGS

# Mavaddat is generating errors

ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/harmonize_alleles.R .
study=Mavaddat_2019_ER_NEG_Breast
# Subset PRS data for each study
awk -v study=$study '$6==study' all_cancer.txt > prs_out/all_cancer.txt_${study}

# Check for duplicate variants based on chr:pos
awk 'a[$1":"$2]++' prs_out/all_cancer.txt_$study | wc -l
# Look for directly matching variants in the WGS data; Remove variants from Mavaddat 2019 that are duplicated and are rare in gnomAD
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/all_cancer.txt_$study plink_data/merged.dat.bim \
| awk '($4==$6 || $4==$7) && ($5==$6 || $5==$7)'| egrep -v 'chr1:51001424:C:CT|chr1:172359627:TA:T|chr2:39472369:CT:C|chr3:49672479:CT:C|chr5:44508162:GT:G|chr5:56366713:CT:C|chr6:20537614:C:CA|chr6:81553832:A:AAT|chr7:140243902:C:CT|chr8:17930101:C:CT|chr10:22188847:A:ACC|chr10:93532430:C:CAA|chr17:45134972:CT:C|chr22:40508703:C:CT' > prs_out/all_cancer.txt_${study}_direct_match
# No direct match
awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $1, $2, $4, $5, $6, a[$1":"$4]}' prs_out/all_cancer.txt_$study plink_data/merged.dat.bim \
| awk '!(($4==$6 || $4==$7) && ($5==$6 || $5==$7))' | grep -v DEL > prs_out/all_cancer.txt_${study}_no_direct_match
# Exclude those that are already a direct match
awk 'NR==FNR{a[$1":"$3];next}!($1":"$3 in a){print}' prs_out/all_cancer.txt_${study}_direct_match prs_out/all_cancer.txt_${study}_no_direct_match \
> prs_out/all_cancer.txt_${study}_no_direct_match_final

wc -l prs_out/all_cancer.txt_${study}_no_direct_match_final




# check how many in PRS and direct match
wc -l prs_out/all_cancer.txt_${study}_direct_match
grep $study all_cancer.txt| wc -l

# grep -vw chr1:113903258:G:T prs_out/all_cancer.txt_${study}_no_direct_match_final
# Check for duplicate variants
# awk 'a[$1":"$3]++' prs_out/all_cancer.txt_${study}_no_direct_match_final > prs_out/all_cancer.txt_${study}_no_direct_match_duplicates
# Drop one from the duplicate; check for the allele frequency first, and get rid of the rare variant keeping the common one
# egrep -vw 'chr9:108126198:G:A|chr3:30641447:G:C' prs_out/all_cancer.txt_${study}_no_direct_match > prs_out/all_cancer.txt_${study}_no_direct_match_uniq
# grep -vw chr9:108126198:G:A prs_out/all_cancer.txt_${study}_no_direct_match > prs_out/all_cancer.txt_${study}_no_direct_match_uniq
# mv prs_out/all_cancer.txt_${study}_no_direct_match_uniq prs_out/all_cancer.txt_${study}_no_direct_match
# Harmonize no direct match alleles
module load R
Rscript harmonize_alleles.R prs_out/all_cancer.txt_${study}_no_direct_match_final

# Update the alleles
awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/all_cancer.txt_${study}_no_direct_match_final_alleles_harmonized > prs_out/all_cancer.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# awk '($NF==1){ print $3, $6, $7, $8, $9}' prs_out/all_cancer.txt_${study}_no_direct_match_final_alleles_harmonized | grep -vw chr9:108126198:G:A > prs_out/all_cancer.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# grep -v chr1:145902073:G:GA prs_out/all_cancer.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt > prs_out/t1
# mv prs_out/t1 prs_out/all_cancer.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt
# Extract study-specific variants
awk '{print $2}' prs_out/all_cancer.txt_${study}_direct_match > prs_out/all_cancer.txt_${study}_direct_match_to_extract.txt
module load plink/1.90b
# If all are direct match, replace prs_out/${study}_direct_match with prs_out/$study and skip next two lines after this line
plink --bfile plink_data/merged.dat --extract prs_out/all_cancer.txt_${study}_direct_match_to_extract.txt --make-bed --out prs_out/${study}_direct_match 
plink --bfile plink_data/merged.dat --extract prs_out/all_cancer.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --update-alleles prs_out/all_cancer.txt_${study}_no_direct_match_alleles_harmonized_update_alleles.txt --make-bed --out prs_out/${study}_harmonized
plink --bfile prs_out/${study}_direct_match --bmerge prs_out/${study}_harmonized --make-bed --out prs_out/$study
# Update variant names
awk '{print $2, $1":"$4}' prs_out/${study}.bim > prs_out/${study}_update_variantnames
plink --bfile prs_out/$study --update-name prs_out/${study}_update_variantnames --make-bed --out prs_out/${study}_varname_updated
# Create a score file
awk '{print $1":"$2, $4, $5}' prs_out/all_cancer.txt_$study > prs_out/${study}.prsweight
# Calculate PRS
plink --bfile prs_out/${study}_varname_updated --score prs_out/${study}.prsweight --out prs_out/${study}_prs


# ## For errors with Mavaddat variants, first check the frequency
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data
# plink --bfile merged.dat --freq --out merged.data.freq
# awk '{print $1":"$3}' prs_out/all_cancer.txt_${study}_direct_match| uniq -c|sed 's/^[ \t]*//;s/[ \t]*$//'| grep ^2
# 2 1:51001424
# 2 1:172359627
# 2 2:39472369
# 2 3:49672479
# 2 5:44508162
# 2 5:56366713
# 2 6:20537614
# 2 6:81553832
# 2 7:140243902
# 2 8:17930101
# 2 10:22188847
# 2 10:93532430
# 2 17:45134972
# 2 22:40508703

# # Now check them individually and remove those that are not true
# grep 22:40508703 merged.data.freq.frq
# grep 22:40508703 ../prs_out/all_cancer.txt_${study}_direct_match
#####################
## CCSS _ Original ##
#####################
module load bcftools/1.9
module load plink/1.90b

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38

ln -s ../chr*.dose.vcf.gz* .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/all_cancer.txt .

# awk 'BEGIN { FS="\t"; OFS="\t" } { $1="chr"$1; $2=$2-1 "\t" $2 } 1' all_cancer.txt| sed '1d' > all_cancer1.bed
# awk 'BEGIN { FS="\t"; OFS="\t" } { $1="chr"$1 } 1' all_cancer.txt| sed '1d' > all_cancer1.bed
# awk 'BEGIN { FS="\t"; OFS="\t" } { $1="chr"$1; $2=$2-1 "\t" $2 } 1' all_cancer.txt| cut -d$'\t' -f1-3,1-3| sed '1d' > all_cancer1.bed

awk '{print "chr"$1"\t"$2-1"\t"$2"\t""chr"$1":"$2-1"-"$2 }' all_cancer.txt| sed '1d' > all_cancer1.bed


# Instead of lifting over the VCF, I will just lift over the PRS SNPs.
/home/aneupane/liftover/liftOver  all_cancer1.bed /home/aneupane/liftover/hg38ToHg19.over.chain all_cancer_GrCh37.bed all_cancer_unmapped.bed

mkdir -p /attr_fraction/plink_data


IFS=$'\n' # set IFS
for line in $(cat all_cancer_GrCh37.bed| sed '1d'); do
VAR="$(echo ${line}| awk '{print $1":"$2"-"$3}'|sed 's/chr//g')"
CHR="$(echo $VAR |awk -F':' '{print $1}')"
echo "Doing ${VAR}"
bcftools view /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/chr${CHR}.dose.vcf.gz ${VAR} > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/plink_data/PRS_${VAR}.vcf.gz
# bcftools query -f '%CHROM:%POS:%REF:%ALT[\n%SAMPLE=%GT]\n' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/yutaka_request_09_27_2022/plink_data/yutaka_${VAR}.vcf.gz >> /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/yutaka_request_09_27_2022/plink_data/yutaka_${VAR}.geno.txt
# annotate VCF
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/plink_data/PRS_${VAR}.vcf.gz -Oz -o /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/plink_data/PRS_${VAR}_edited.vcf.gz
plink --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/plink_data/PRS_${VAR}_edited.vcf.gz --double-id --vcf-half-call m --keep-allele-order --make-bed --out /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/plink_data/PRS_${VAR}_edited 2>&1 | tee -a extract_plink_all.log
done



bcftools view /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/chr1.dose.vcf.gz 1:121287993-121287994| less -S
1:10566214-10566215
1:110198128-110198129

