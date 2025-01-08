# Note:: In /biallelic, I did not apply minGQ and min-meanDP filters, but in biallelic2, I did apply those filters.
######################################################################################

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC

hpcf_interactive -n 12 -R "rusage[mem=8000]" -q standard
# ## Convert to biallelic
# bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa MERGED_sorted_sjlife_1_2_zhaoming_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_zhaoming_v2.vcf.gz
# bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' MERGED_biallelic_sorted_sjlife_1_2_zhaoming_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz
# bcftools index -f -t --threads 4 MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz

ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/*SJLIFE_CCSS_WES_101724.GATK4180.hg38.vcf.gz* .

## Incase we need to rename
awk '{cmd="echo "$0" | sed -e '\''s/.*CCSS-//g; s/^0\\+\\([^0]\\)/\\1/g'\''"; cmd | getline result; close(cmd); print $0"\t"result}' Survivor_WES.samplelist.ccss  > ./biallelic2/rename_ccss.txt
cat rename_ccss.txt rename_sjlife.txt > sample_mapping.txt

## <biallelic_renaming_split.sh>
for i in {1..22} X Y; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
unset VCF; \
export THREADS=4; \
export VCF="${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38.vcf.gz"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/"; \
bsub \
        -P "${CHR}_biallelic" \
        -q "standard" \
        -J "${CHR}_biallelic" \
        -o "${WORKDIR}/biallelic2/logs/${VCF%.vcf*}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/biallelic_renaming_split.sh"; \
done;


zcat chr1.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.vcf.gz | head -10000 | grep "^#CHROM" | awk '{for (i=10; i<=NF; i++) print $i}' > VCFsample_names.txt

## Analysis version 2 (QCed for GQ, DP and VQSR )
## Applying GQ filter
# <biallelic_renaming_split2.sh>
#!/usr/bin/bash
module load bcftools
module load vcftools/0.1.16
module load tabix
vcftools --gzvcf "${WORKDIR}/${VCF}" --chr ${CHR} --keep-filtered PASS --minGQ 20 --min-meanDP 10 --recode --stdout | bgzip > "${WORKDIR}/biallelic2/$(basename ${VCF} .vcf.gz)_tmp0.vcf.gz"
# ## Normalize and make biallelic
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa "${WORKDIR}/biallelic2/$(basename ${VCF} .vcf.gz)_tmp0.vcf.gz" -Oz -o "${WORKDIR}/biallelic2/$(basename ${VCF} .vcf.gz)_tmp.vcf.gz"
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' "${WORKDIR}/biallelic2/$(basename ${VCF} .vcf.gz)_tmp.vcf.gz" -Oz -o "${WORKDIR}/biallelic2/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
bcftools index -f -t --threads 4 "${WORKDIR}/biallelic2/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"


# # ## Extract three cohorts
# bcftools view -O z -o \
# "${WORKDIR}/biallelic/ccss/$(basename ${VCF} .vcf.gz)_biallelic_ccss.vcf.gz" \
# -S ${WORKDIR}/biallelic/extract_CCSS.samples.txt \
# "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
# bcftools index -f -t --threads 4 "${WORKDIR}/biallelic/ccss/$(basename ${VCF} .vcf.gz)_biallelic_ccss.vcf.gz"

# bcftools view -O z -o \
# "${WORKDIR}/biallelic/sjlife/$(basename ${VCF} .vcf.gz)_biallelic_sjlife.vcf.gz" \
# -S ${WORKDIR}/biallelic/extract_SJLIFE_survivor.txt \
# "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
# bcftools index -f -t --threads 4 "${WORKDIR}/biallelic/sjlife/$(basename ${VCF} .vcf.gz)_biallelic_sjlife.vcf.gz"

# bcftools view -O z -o \
# "${WORKDIR}/biallelic/sjlife_control/$(basename ${VCF} .vcf.gz)_biallelic_sjlife_control.vcf.gz" \
# -S ${WORKDIR}/biallelic/extract_SJLIFE_survivor_control.txt \
# "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
# bcftools index -f -t --threads 4 "${WORKDIR}/biallelic/sjlife_control/$(basename ${VCF} .vcf.gz)_biallelic_sjlife_control.vcf.gz"


## # Rename samples (if we need to rename)
## bcftools reheader -s sample_mapping.txt "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz" -o "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic_renamed.vcf.gz"


## Remove LCR region
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/plink_all
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4271055/
# wget https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs38.bed.gz
ln -s ../biallelic/plink_all/LCR-hs38.bed .

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/
for i in {1..22} X Y; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
unset VCF; \
export THREADS=4; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/"; \
bsub \
        -P "${CHR}_plk" \
        -J "${CHR}_plk" \
        -o "${WORKDIR}/biallelic2/plink_all/logs/${VCF%.vcf*}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "./plinkQC.sh"; \
done;


# <plinkQC.sh>
#!/usr/bin/bash
## Variant level QC
# Filter based on call rate (<90%)
module load plink/1.90b
plink --vcf "${WORKDIR}/biallelic2/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.vcf.gz" --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic"
plink --bfile "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic" --geno 0.10 --keep-allele-order --make-bed --out "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1"
## Chromosome Y
# Before main variant filters, 8065 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is -nan.
# Error: All variants excluded due to missing genotype data (--geno).


# Filter based on HWE p-value (<1x10^-15)
plink --bfile "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1" --hwe 1e-15 --keep-allele-order --make-bed --out "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15"
# Exclude variants in low-complexity regions (use BED file)
# bedtools subtract -a filtered_hwe_data.bim -b low_complexity_regions.bed > final_filtered_data.bim
# plink --bfile final_filtered_data --make-bed --out final_filtered_data
plink --bfile "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15" --exclude range /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/LCR-hs38.bed --keep-allele-order --make-bed --out "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed"
# Calculate MAC
module load plink/2.0
plink2 --bfile "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed" --mac 1 --keep-allele-order --make-bed --out "${WORKDIR}/biallelic2/plink_all/${CHR}.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1"


## Merge original plink to get variant counts
ls chr*.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.bim | sed 's/.bim//'| sort -V > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic


plink --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic --update-ids /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt --make-bed --keep-allele-order --out chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated



## This needs to be updated
# <plinkQC.sh>
#!/usr/bin/bash
## Variant level QC
# Filter based on call rate (<90%)
Biallelic: 2011227
GQ, DP and VQSR: 1523098
90% call rate: 1356453
HWE: 1350201
LCR removed: 1342006
MAC greater or equal 1: 1340181



## sample level QC
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/plink_all
ls chr*.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.bim | sed 's/.bim//'| sort -V  > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15


# plink --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 --indep-pairwise 100 25 0.2 --keep-allele-order --maf 0.05 --make-bed --out chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune


# plink \
#   --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 \
#   --indep-pairwise 100 25 0.2 \
#   --keep-allele-order \
#   --maf 0.05 \
#   --out geno.0.1.hwe.1e-15_indep_pairwise

plink \
  --bfile chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 \
  --exclude range /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/high-LD-regions-GRCh38.txt \
  --indep-pairwise 100 25 0.2 \
  --keep-allele-order \
  --maf 0.05 \
  --out geno.0.1.hwe.1e-15_indep_pairwise

plink --bfile chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 \
--keep-allele-order \
--extract geno.0.1.hwe.1e-15_indep_pairwise.prune.in \
--make-bed \
--out chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune

plink --bfile chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune --het --out chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune_heterozygosity
plink --bfile chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 --missing --out chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_missing


## Merge plink
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all
ls chr*.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.bim | sed 's/.bim//'| sort -V  > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1

## Update IDs
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --update-ids /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt --make-bed --keep-allele-order --out chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated



#######################################################
## Extract Survivorr, control and CCSS exp here from ##
#######################################################
## Survivor; ../extract_SJLIFE_survivor_iid_fid.txt
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --keep-allele-order --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_SJLIFE_survivor_iid_fid.txt --make-bed --out Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1
plink --bfile  Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --update-ids /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt --make-bed --keep-allele-order --out Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated
## SJLIFE als has 35 duplicates, we can remove them based on low call rates:
plink --bfile chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated --missing --out chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated_missing
plink --bfile chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated --remove duplicate_samples_to_remove.txt --keep-allele-order --make-bed --out chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated_unique

## CCSS exp; ../extract_CCSS.samples_iid_fid.txt
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --keep-allele-order --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_CCSS.samples_iid_fid.txt --make-bed --out CCSS_exp/chr.ALL.CCSS_exp_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1
plink --bfile CCSS_exp/chr.ALL.CCSS_exp_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --update-ids /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt --make-bed --keep-allele-order --out CCSS_exp/chr.ALL.CCSS_exp_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated
## Control; ../extract_SJLIFE_survivor_control_iid_fid.txt
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --keep-allele-order --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_SJLIFE_survivor_control_iid_fid.txt --make-bed --out Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1
plink --bfile Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --update-ids /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt --make-bed --keep-allele-order --out Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated


## We can still use the same annotation done on cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation
## Install LOFTEE
https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache
https://wiki.stjude.org/display/CAB/Ensembl+-+VEP
https://wiki.stjude.org/pages/viewpage.action?pageId=40110325

## Revel: https://sites.google.com/site/revelgenomics/downloads; details: https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#revel
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv
## for GRCh37:
## tabix -f -s 1 -b 2 -e 2 new_tabbed_revel.tsv.gz
## for GRCh38:
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz






cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff/revel_loftee
ln -s ../new_chr*.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf .







# ## Those that failed to generate FIELDS
# for i in 2 3 4 5 6 7 8 9 10 11 14 15 16 17 20 22; do 
# CHR="chr${i}"; 
# echo "DOing ${CHR}"
# ANNOTATED="new_${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp_clinvar_12_10_2023.vcf"
# cat ${ANNOTATED}|${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF" "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_aaalt" "dbNSFP_aapos" "dbNSFP_aaref" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_raw" "AF" > ${WORKDIR}/"$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf-FIELDS-simple.txt"
# done





# ## Basic
# # --vcf; --tab
# vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
# --fasta /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
# --offline --everything --assembly GRCh38 \
# --dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ \
# -i new_chr18.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp.vcf \
# --tab \
# -o test.tsv \
# --plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:false
# # vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache --offline --everything --assembly GRCh38 --dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ --fasta /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa --force_overwrite -i new_chr8.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf --tab -o VEP.GRCh38.with.gnomAD.tsv --plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:false
# # vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache --offline --everything --vcf --assembly GRCh38 --dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ --fasta /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa --force_overwrite -i new_chr8.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf --tab -o VEP.GRCh38.with.gnomAD.tsv --plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/human_ancestor.fa.gz

# vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
# --fasta /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
# --offline --everything --fields ALL --assembly GRCh38 \
# -i new_chr18.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp.vcf \
# --tab \
# -o test.tsv 


## Calculate mean DP
## <calculate_mean_coverage.sh>


#!/bin/bash

# Directory containing VCF files for each chromosome
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/

#!/bin/bash
module load bcftools
# Directory containing your VCF files
vcf_directory="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/"
output_file="mean_coverage.txt"

# Initialize an array to store mean coverages
mean_coverages=()

# Loop through chromosomes 1 to 22
# Number of samples
total_samples=8073

for chromosome in {1..22}; do
    vcf_file="${vcf_directory}chr${chromosome}.Survivor_WES.GATK4180.hg38_biallelic.vcf.gz"

    if [ -f "$vcf_file" ]; then
        # Extract DP values and calculate the mean coverage for the current chromosome
        mean_coverage_tmp=$(bcftools query -f '%INFO/DP\n' "$vcf_file" | awk '{sum+=$1}END{print sum/NR}')
        # Divide the mean coverage by the total number of samples
        mean_coverage=$(echo "scale=2; $mean_coverage_tmp / $total_samples" | bc)
        echo "Chromosome $chromosome: Mean Coverage = $mean_coverage"
        mean_coverages+=("$mean_coverage")
    fi
done

# Calculate the mean coverage across all chromosomes
total_mean_coverage=0
for coverage in "${mean_coverages[@]}"; do
    total_mean_coverage=$(echo "$total_mean_coverage + $coverage" | bc)
done
num_chromosomes=${#mean_coverages[@]}
overall_mean_coverage=$(echo "scale=2; $total_mean_coverage / $num_chromosomes" | bc)

# Save the results to an output file
echo "Mean Coverage Across All Chromosomes: $overall_mean_coverage" > "$output_file"
echo "Individual Mean Coverages:" >> "$output_file"
for chromosome in {1..22}; do
    vcf_file="${vcf_directory}chr${chromosome}.Survivor_WES.GATK4180.hg38_biallelic.vcf.gz"
    if [ -f "$vcf_file" ]; then
        mean_coverage_tmp=$(bcftools query -f '%INFO/DP\n' "$vcf_file" | awk '{sum+=$1}END{print sum/NR}')
        # Divide the mean coverage by the total number of samples
        mean_coverage=$(echo "scale=2; $mean_coverage_tmp / $total_samples" | bc)
        echo "Chromosome $chromosome: Mean Coverage = $mean_coverage" >> "$output_file"
    fi
done

echo "Results saved to $output_file"


head MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.PASS.decomposed.bim |cut -d$'\t' -f2|awk -F':' '{split($3 a $4, alleles, ""); asort(alleles); printf("%s:%s:%s:%s\n", $1, $2, alleles[1], alleles[2])}'

# ## Calculate mean DP
# ## <calculate_mean_coverage2.sh>

# #!/bin/bash
# module load bcftools
# # Directory containing your VCF files
# vcf_directory="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/"
# output_file="mean_coverage2.txt"

# # Initialize variables to store total coverage and total variants
# total_coverage=0
# total_variants=0

# # Function to format a number without scientific notation
# format_number() {
#     local number="$1"
#     printf "%.0f" "$number"
# }

# # Loop through chromosomes 1 to 22
# for chromosome in {1..22}; do
#     vcf_file="${vcf_directory}chr${chromosome}.Survivor_WES.GATK4180.hg38_biallelic.vcf.gz"

#     if [ -f "$vcf_file" ]; then
#         # Use bcftools to get the total variant count and total DP (coverage) for the chromosome
#         chromosome_variants=$(bcftools view -H -c1 "$vcf_file" | wc -l)
#         chromosome_coverage=$(bcftools query -f '%INFO/DP\n' "$vcf_file" | awk -v total_samples="8073" '{sum += $1 / total_samples} END {print sum}')
        
#         # Add the chromosome's coverage and variants to the totals
#         total_coverage=$(echo "$total_coverage + $chromosome_coverage" | bc)
#         total_variants=$((total_variants + chromosome_variants))

#         # Print and save the results for each chromosome
#         echo "Chromosome ${chromosome}:"
#         formatted_chromosome_coverage=$(format_number "${chromosome_coverage}")
#         formatted_chromosome_variants=$(format_number "${chromosome_variants}")
#         echo "Total Coverage: ${formatted_chromosome_coverage}"
#         echo "Total Variants: ${formatted_chromosome_variants}"

#         # Save the results in a separate file for each chromosome
#         echo "Total Coverage in chr ${chromosome}: ${formatted_chromosome_coverage}" >> "$output_file"
#         echo "Total Variants in chr ${chromosome}: ${formatted_chromosome_variants}" >> "$output_file"
#     fi
# done

# # Calculate final mean coverage
# if [ "$total_variants" -gt 0 ]; then
#     final_mean_coverage=$(echo "scale=2; $total_coverage / $total_variants" | bc)
#     formatted_final_mean_coverage=$(format_number "${final_mean_coverage}")
#     echo "Final mean Coverage: ${formatted_final_mean_coverage}"
# else
#     total_coverage_per_variant=0
# fi

# # Save the total results in the final output file
# formatted_total_coverage=$(format_number "${total_coverage}")
# formatted_total_variants=$(format_number "${total_variants}")
# echo "Total Coverage Across All Chromosomes: ${formatted_total_coverage}" >> "$output_file"
# echo "Total Variants Across All Chromosomes: ${formatted_total_variants}" >> "$output_file"
# echo "Final mean Coverage: ${formatted_final_mean_coverage}" >> "$output_file"
# echo "Results saved to ${output_file}"








