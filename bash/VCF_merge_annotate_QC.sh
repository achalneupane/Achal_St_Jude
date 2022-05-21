#!/bin/bash
### ##################
## Date: 04/22/2022 ##
## Achal Neupane    ##
######################
#################################
## Phenotype data manipulation ##
#################################
# Wrote an R code : phenotype_cleaning_st_jude_life.r

##########################################################
## Script to merge VCF files from disjoint sample lists ##
##########################################################
# ln -s ../sjlife_2/SJLIFE2.GATKv3.4.VQSR.sjlid_chr22.PASS.decomposed.vcf.gz* .
# ln -s ../sjlife_2/SJLIFE2.GATKv3.4.VQSR.sjlid_chr20.PASS.decomposed.vcf.gz* .
# ln -s ../sjlife_1/SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr22.PASS.decomposed.sjlid.qced.vcf.gz* .
# ln -s ../sjlife_1/SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr20.PASS.decomposed.sjlid.qced.vcf.gz* .

ln -s ../sjlife_2/SJLIFE2.GATKv3.4.VQSR.sjlid_chr*.PASS.decomposed.vcf.gz* .
ln -s ../sjlife_1/SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr*.PASS.decomposed.sjlid.qced.vcf.gz* .
(module avail) |& grep -i bcftools| sort -V

for CHR in {1..22}; do \
	unset VCF1; unset VCF2; unset MERGED; \
	echo "Doing chr${CHR}"; \
	export VCF1="SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${CHR}.PASS.decomposed.sjlid.qced.vcf.gz"; \
	export VCF2="SJLIFE2.GATKv3.4.VQSR.sjlid_chr${CHR}.PASS.decomposed.vcf.gz"; \
	echo -e "**\nMerging ${VCF1} \nand \n${VCF2}\n**"; \
	export MEM=6; \
	export THREADS=4; \
	export MERGED="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed"; \
	export OUT_DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/"; \
	bsub \
	-P "chr${CHR}_merge" \
	-J "chr${CHR}_merge" \
	-o "${OUT_DIR}/logs/${MERGED}_s00VCFmerge.%J" \
	-n ${THREADS} \
	-R "rusage[mem=8192]" \
	"./entrypoint_merge_vcf.sh"; \
done; 


#############################
## VCF to plink conversion ##
#############################
for CHR in {1..22}; do \
	unset VCF;unset PLINK_FILE; \
	echo "Doing chr${CHR}"; \
	export VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz"; \
	export PLINK_FILE="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${CHR}.PASS.decomposed"; \
	export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/"
	echo -e "**\nConverting ${VCF} \n to Plink \n${PLINK_FILE}\n**"; \
	export MEM=6; \
	export THREADS=4; \
	bsub \
	-P "chr${CHR}_plink" \
	-J "chr${CHR}_plink" \
	-o "${WORKDIR}/logs/${PLINK_FILE}_s01VCFplink.%J" \
	-n ${THREADS} \
	-R "rusage[mem=8192]" \
	"./entrypoint_VCF_to_Plink_with_basic_QC.sh"; \
done; 

###############
## Admixture ##
###############
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/1kGP/1000genomes_merged.* .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise.* .

# 1000G
# first fix bim file of 1000G
cut -f2 1000genomes_merged.bim > 1KGsnps
cut -f2 SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise.bim > SJLIFEsnps
# Find common SNPs between 1KG and SJLIFE
cat 1KGsnps SJLIFEsnps|sort |uniq -c |sed -n -e 's/^ *2 \(.*\)/\1/p' > commonSNPs_in_1KG_SJLIFE


## Extract common SNPs from 1KG
plink1.9 --memory 300000 --threads 24 --bfile ${BFILE} --extract commonSNPs_in_1KG_NHW --geno 0.02 --maf 0.02 --keep-allele-order --make-bed --out ${BFILE}_with_1KGsnps

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
####################
## VCF annotation ##
####################
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation
ln -s ../*.vcf.gz .

## Find column number of MetaSVM_pred
# C=1
# for i in $(zcat /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz| head -n 1) ; do
#     echo $i
#     if [ $i == "MetaSVM_pred" ] ; then
#         break ;
#     else
#         echo $C
#         C=$(( $C + 1 ))
#     fi
# done

# bgzip -c test.vcf > test.vcf.gz
# bcftools index -t --threads 4 test.vcf
# bgzip -c> ${WXS}_with_no_chr.vcf.gz && tabix -s1 -b2 -e2 ${VCF}_with_no_chr.vcf.gz
#####################################################################
## annotate dbSNP VCF to add chr string infront of the chromosomes ##
#####################################################################
## Downloaded latest dbSNP on 05/02/2022 using  wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.39.gz; wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.39.gz.tbi
#### Downloaded latest dbSNP on 05/02/2022 using  wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
# First check the chromosomes in VCF
zcat GCF_000001405.39.gz|grep -v "^#"|cut -d$'\t' -f1| uniq

cat << EoF > chrchange.txt
NC_000001.11 chr1
NC_000002.12 chr2
NC_000003.12 chr3
NC_000004.12 chr4
NC_000005.10 chr5
NC_000006.12 chr6
NC_000007.14 chr7
NC_000008.11 chr8
NC_000009.12 chr9
NC_000010.11 chr10
NC_000011.10 chr11
NC_000012.12 chr12
NC_000013.11 chr13
NC_000014.9	chr14
NC_000015.10 chr15
NC_000016.10 chr16
NC_000017.11 chr17
NC_000018.10 chr18
NC_000019.10 chr19
NC_000020.11 chr20
NC_000021.9 chr21
NC_000022.11 chr22
NC_000023.11 chrX
NC_000024.10 chrY
NC_012920.1 chrMT
EoF

bcftools annotate --rename-chrs chrchange.txt GCF_000001405.39.gz -Oz -o dbSNP_155.GCF_000001405.39.gz

bcftools index -t --threads 4 dbSNP_155.GCF_000001405.39.gz

# ## BFILE with 1502 individuals
# BFILE="FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC2"
# ## Plink to VCF file conversion
# VCF="FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC2_VCF"
# plink1.9 --bfile ${BFILE} --keep-allele-order --allow-no-sex --recode vcf-iid --output-chr MT --out ${VCF}
# Total genotyping rate is 0.997404.
# 2901993 variants and 1502 people pass filters and QC.
## Download dbnsfp (05/03/2022)
# wget https://snpeff.blob.core.windows.net/databases/dbs/GRCh38/dbNSFP_4.1a/dbNSFP4.1a.txt.gz; wget https://snpeff.blob.core.windows.net/databases/dbs/GRCh38/dbNSFP_4.1a/dbNSFP4.1a.txt.gz.tbi
## from dbnsfp "If you used our ensemble scores MetaSVM and MetaLR, which are based on 10 component scores (SIFT, PolyPhen-2 HDIV, PolyPhen-2 HVAR, GERP++, MutationTaster, Mutation Assessor, FATHMM, LRT, SiPhy, PhyloP) and the maximum frequency observed in the 1000 genomes populations. Please cite:
## 1. Dong C, Wei P, Jian X, Gibbs R, Boerwinkle E, Wang K* and Liu X*. (2015) Comparison and integration of deleteriousness prediction methods for nonsynonymous SNVs in whole exome sequencing studies. Human Molecular Genetics 24(8):2125-2137. *corresponding authors [PDF]"
## Definition of terms (eg., MetaSVM_pred) https://github.com/achalneupane/st_jude_code/blob/master/references/dbnsfpv3_preprint.pdf

## List of databases 
# java -jar snpEff.jar databases| grep something


# # Exac GrCh38
# downloaded Exac database on 05/02/2022; wget http://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/ExAC.0.3.GRCh38.vcf.gz

## Download clinvar
# Date: 05/04/2022; wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz; wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi






# /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/table_annovar.pl test.intervar.vcf \
# /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/humandb \
# -buildver hg38 \
# -out myanno -remove \
# -protocol refGene,exac03,exac03nontcga,esp6500siv2_all,gnomad_exome,gnomad_genome,dbnsfp42c,revel,intervar_20180118,dbscsnv11,cosmic70,nci60,clinvar_20220320 \
# -operation g,f,f,f,f,f,f,f,f,f,f,f,f \
# -nastring . -vcfinput 
 
## Generate test VCF
# zcat test.vcf.gz| head -10000 |bgzip -c> test_chr1.vcf.gz &&  tabix -s1 -b2 -e2 test_chr1.vcf.gz






for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.PASS.decomposed.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
	-P "${CHR}_annotate" \
	-J "${CHR}_annotate" \
	-o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
	-n ${THREADS} \
	-R "rusage[mem=8192]" \
	"./entrypoint_VCFannotation.sh"; \
done; 



#####################################
## Annovar and intervar Annotation ##
#####################################
# First check avblist hg38_avdblist.txt for the latest releases
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar
perl annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg38 /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/annovar/humandb/

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar 1000g2015aug humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03nontcga humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar esp6500siv2_all humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_exome humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar revel humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar intervar_20180118 humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbscsnv11 humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic70 humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar nci60 humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20220320 humandb/ 


## Annovar annotation ##
# https://annovar.openbioinformatics.org/en/latest/user-guide/startup/
# NOTE: the -operation argument tells ANNOVAR which operations to use for each of the protocols: g means gene-based, gx means gene-based with cross-reference annotation (from -xref argument), r means region-based and f means filter-based.

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/ANNOVAR_ANNOTATION
ln -s ../../*.vcf.gz .

## Intervar ##
## Add Intervar annotation
# --input_type AVinput or VCF or VCF_m

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.PASS.decomposed.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/ANNOVAR_ANNOTATION/"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
bsub \
-P "${CHR}_ANNOVAR" \
-J "${CHR}_ANNOVAR" \
-o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_ANNOVAR.%J" \
-n ${THREADS} \
-R "rusage[mem=8192]" \
"./entrypoint_ANNOVAR_annotation.sh"; \
done;



##################
## PCA analysis ##
##################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA
ln -s ../final.bed .; ln -s ../final.bim .; ln -s ../final.fam .
## Population file
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/1kGP/integrated_call_samples_v3.20130502.ALL.panel .


module load plink/1.90b
BFILE="final"
plink --memory 400000 --threads 24 --bfile ${BFILE} --allow-no-sex --cluster --geno 0.01 --genome --hwe 0.001 --ld-window-r2 0.2 --maf 0.02 --mds-plot 4 --min 0.2 --nonfounders --pca header --out ${BFILE}-PCAS
# 1031903 MB RAM detected; reserving 400000 MB for main workspace.
# 188445 variants loaded from .bim file.
# 6985 people (0 males, 0 females, 6985 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final-PCAS.nosex .
# Using up to 24 threads (change this with --threads).
# Warning: This run includes BLAS/LAPACK linear algebra operations which
# currently disregard the --threads limit.  If this is problematic, you may want
# to recompile against single-threaded BLAS/LAPACK.
# Before main variant filters, 6985 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.999482.
# 2 variants removed due to missing genotype data (--geno).
# --hwe: 48310 variants removed due to Hardy-Weinberg exact test.
# 0 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 140133 variants and 6985 people pass filters and QC.
# Note: No phenotypes present.
# Relationship matrix calculation complete.
# --pca: Results saved to final-PCAS.eigenval and final-PCAS.eigenvec .
# IBD calculations complete.
# Finished writing final-PCAS.genome .
# Clustering... done.
# Cluster solution written to final-PCAS.cluster1 , final-PCAS.cluster2 , and
# final-PCAS.cluster3 .
# Performing multidimensional scaling analysis (SVD algorithm, 4
# dimensions)... done.
# MDS solution written to final-PCAS.mds .

@@ Now I am using an R script I wrote to plot and select samples: PCA_analysis_SJLIFE.r

#########################
## Cleaned PCA round 1 ##
#########################
## Extract samples for PCA with three ethnicities only
BFILE="final"
plink --memory 300000 --threads 24 --bfile ${BFILE} --keep samples.to.exclude.round2.pca.txt --make-bed --keep-allele-order --out ${BFILE}_cleaned1
# 188445 variants loaded from .bim file.
# 6985 people (0 males, 0 females, 6985 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_cleaned1.nosex .
# --keep: 5754 people remaining.
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 5754 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate in remaining samples is 0.999419.
# 188445 variants and 5754 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to final_cleaned1.bed + final_cleaned1.bim + final_cleaned1.fam ...

BFILE="final_cleaned1"
plink --memory 400000 --threads 24 --bfile ${BFILE} --allow-no-sex --cluster --geno 0.01 --genome --hwe 0.001 --ld-window-r2 0.2 --maf 0.02 --mds-plot 4 --min 0.2 --nonfounders --pca header --out ${BFILE}-PCAS
# 188445 variants loaded from .bim file.
# 5754 people (0 males, 0 females, 5754 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_cleaned1-PCAS.nosex .
# Using up to 24 threads (change this with --threads).
# Warning: This run includes BLAS/LAPACK linear algebra operations which
# currently disregard the --threads limit.  If this is problematic, you may want
# to recompile against single-threaded BLAS/LAPACK.
# Before main variant filters, 5754 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.999419.
# 2 variants removed due to missing genotype data (--geno).
# --hwe: 35259 variants removed due to Hardy-Weinberg exact test.
# 0 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 153184 variants and 5754 people pass filters and QC.
# Note: No phenotypes present.
# Relationship matrix calculation complete.
# --pca: Results saved to final_cleaned1-PCAS.eigenval and
# final_cleaned1-PCAS.eigenvec .


#########################
## Cleaned PCA round 2 ##
#########################

## Round 2 for per ethnicity from PCA
## Extract samples for second round of PCA
## EUR
BFILE="final"
plink --memory 300000 --threads 24 --bfile ${BFILE} --keep SJLIFE_EUR_Per_PCA.txt --make-bed --keep-allele-order --out ${BFILE}_EUR
# 188445 variants loaded from .bim file.
# 6985 people (0 males, 0 females, 6985 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_EUR.nosex .
# --keep: 3439 people remaining.
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 3439 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate in remaining samples is 0.999186.
# 188445 variants and 3439 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to final_EUR.bed + final_EUR.bim + final_EUR.fam ... done.


BFILE="final"
plink --memory 300000 --threads 24 --bfile ${BFILE} --keep SJLIFE_AFR_Per_PCA.txt --make-bed --keep-allele-order --out ${BFILE}_AFR
# 188445 variants loaded from .bim file.
# 6985 people (0 males, 0 females, 6985 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_AFR.nosex .
# --keep: 645 people remaining.
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 645 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate in remaining samples is 0.999168.
# 188445 variants and 645 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to final_AFR.bed + final_AFR.bim + final_AFR.fam ... done.



BFILE="final"
plink --memory 300000 --threads 24 --bfile ${BFILE} --keep SJLIFE_EAS_Per_PCA.txt --make-bed --keep-allele-order --out ${BFILE}_EAS
# 188445 variants loaded from .bim file.
# 6985 people (0 males, 0 females, 6985 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_EAS.nosex .
# --keep: 16 people remaining.
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 16 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate in remaining samples is 0.99966.
# 188445 variants and 16 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to final_EAS.bed + final_EAS.bim + final_EAS.fam ... done.



BFILE="final_EUR"
plink --memory 400000 --threads 24 --bfile ${BFILE} --allow-no-sex --cluster --geno 0.05 --genome --hwe 1e-06 --ld-window-r2 0.2 --maf 0.05 --mds-plot 4 --min 0.2 --nonfounders --pca header --out ${BFILE}-PCAS_EUR
# 188445 variants loaded from .bim file.
# 3439 people (0 males, 0 females, 3439 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_EUR-PCAS_EUR.nosex .
# Using up to 24 threads (change this with --threads).
# Warning: This run includes BLAS/LAPACK linear algebra operations which
# currently disregard the --threads limit.  If this is problematic, you may want
# to recompile against single-threaded BLAS/LAPACK.
# Before main variant filters, 3439 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.999186.
# 0 variants removed due to missing genotype data (--geno).
# --hwe: 30 variants removed due to Hardy-Weinberg exact test.
# 7742 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 180673 variants and 3439 people pass filters and QC.
# Note: No phenotypes present.
# Relationship matrix calculation complete.
# --pca: Results saved to final_EUR-PCAS_EUR.eigenval and
# final_EUR-PCAS_EUR.eigenvec .


BFILE="final_AFR"
plink --memory 400000 --threads 24 --bfile ${BFILE} --allow-no-sex --cluster --geno 0.05 --genome --hwe 1e-06 --ld-window-r2 0.2 --maf 0.05 --mds-plot 4 --min 0.2 --nonfounders --pca header --out ${BFILE}-PCAS_AFR
# 188445 variants loaded from .bim file.
# 645 people (0 males, 0 females, 645 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_AFR-PCAS_AFR.nosex .
# Using up to 24 threads (change this with --threads).
# Warning: This run includes BLAS/LAPACK linear algebra operations which
# currently disregard the --threads limit.  If this is problematic, you may want
# to recompile against single-threaded BLAS/LAPACK.
# Before main variant filters, 645 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.999168.
# 20 variants removed due to missing genotype data (--geno).
# --hwe: 63 variants removed due to Hardy-Weinberg exact test.
# 53500 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 134862 variants and 645 people pass filters and QC.
# Note: No phenotypes present.
# Relationship matrix calculation complete.
# --pca: Results saved to final_AFR-PCAS_AFR.eigenval and
# final_AFR-PCAS_AFR.eigenvec .
# IBD calculations complete.
# Finished writing final_AFR-PCAS_AFR.genome .
# Clustering... done.
# Cluster solution written to final_AFR-PCAS_AFR.cluster1 ,
# final_AFR-PCAS_AFR.cluster2 , and final_AFR-PCAS_AFR.cluster3 .
# Performing multidimensional scaling analysis (SVD algorithm, 4
# dimensions)... done.


BFILE="final_EAS"
plink --memory 400000 --threads 24 --bfile ${BFILE} --allow-no-sex --cluster --geno 0.05 --genome --hwe 1e-06 --ld-window-r2 0.2 --maf 0.05 --mds-plot 4 --min 0.2 --nonfounders --pca header --out ${BFILE}-PCAS_EAS
# 188445 variants loaded from .bim file.
# 16 people (0 males, 0 females, 16 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_EAS-PCAS_EAS.nosex .
# Using up to 24 threads (change this with --threads).
# Warning: This run includes BLAS/LAPACK linear algebra operations which
# currently disregard the --threads limit.  If this is problematic, you may want
# to recompile against single-threaded BLAS/LAPACK.
# Before main variant filters, 16 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.99966.
# 953 variants removed due to missing genotype data (--geno).
# --hwe: 0 variants removed due to Hardy-Weinberg exact test.
# 71544 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 115948 variants and 16 people pass filters and QC.
# Note: No phenotypes present.
# Relationship matrix calculation complete.
# Warning: calculating 16 PCs, since there are only 16 samples.
# --pca: Results saved to final_EAS-PCAS_EAS.eigenval and
# final_EAS-PCAS_EAS.eigenvec .
# IBD calculations complete.
# Finished writing final_EAS-PCAS_EAS.genome .
# Clustering... done.
# Cluster solution written to final_EAS-PCAS_EAS.cluster1 ,
# final_EAS-PCAS_EAS.cluster2 , and final_EAS-PCAS_EAS.cluster3 .
# Performing multidimensional scaling analysis (SVD algorithm, 4
# dimensions)... done.
# MDS solution written to final_EAS-PCAS_EAS.mds .

##################################
## Repeat round 2 for Admixture ##
##################################


## Round 2 for per ethnicity from PCA
## Extract samples for second round from Admixture
## EUR
BFILE="final"
plink --memory 300000 --threads 24 --bfile ${BFILE} --keep SJLIFE_EUR_Per_ADMIXTURE.txt --make-bed --keep-allele-order --out ${BFILE}_EUR_ADMIXTURE
# 188445 variants loaded from .bim file.
# 6985 people (0 males, 0 females, 6985 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_EUR_ADMIXTURE.nosex .
# --keep: 3423 people remaining.
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 3423 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate in remaining samples is 0.999185.
# 188445 variants and 3423 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to final_EUR_ADMIXTURE.bed + final_EUR_ADMIXTURE.bim +
# final_EUR_ADMIXTURE.fam ... done.



BFILE="final"
plink --memory 300000 --threads 24 --bfile ${BFILE} --keep SJLIFE_AFR_Per_ADMIXTURE.txt --make-bed --keep-allele-order --out ${BFILE}_AFR_ADMIXTURE
# 188445 variants loaded from .bim file.
# 6985 people (0 males, 0 females, 6985 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_AFR_ADMIXTURE.nosex .
# --keep: 592 people remaining.
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 592 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate in remaining samples is 0.999131.
# 188445 variants and 592 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to final_AFR_ADMIXTURE.bed + final_AFR_ADMIXTURE.bim +




BFILE="final"
plink --memory 300000 --threads 24 --bfile ${BFILE} --keep SJLIFE_EAS_Per_ADMIXTURE.txt --make-bed --keep-allele-order --out ${BFILE}_EAS_ADMIXTURE
# 188445 variants loaded from .bim file.
# 6985 people (0 males, 0 females, 6985 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_EAS_ADMIXTURE.nosex .
# --keep: 21 people remaining.
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 21 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate in remaining samples is 0.999665.
# 188445 variants and 21 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to final_EAS_ADMIXTURE.bed + final_EAS_ADMIXTURE.bim +
# final_EAS_ADMIXTURE.fam ... done.



BFILE="final_EUR_ADMIXTURE"
plink --memory 400000 --threads 24 --bfile ${BFILE} --allow-no-sex --cluster --geno 0.05 --genome --hwe 1e-06 --ld-window-r2 0.2 --maf 0.05 --mds-plot 4 --min 0.2 --nonfounders --pca header --out ${BFILE}-PCAS_EUR
# 188445 variants loaded from .bim file.
# 3423 people (0 males, 0 females, 3423 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_EUR_ADMIXTURE-PCAS_EUR.nosex .
# Using up to 24 threads (change this with --threads).
# Warning: This run includes BLAS/LAPACK linear algebra operations which
# currently disregard the --threads limit.  If this is problematic, you may want
# to recompile against single-threaded BLAS/LAPACK.
# Before main variant filters, 3423 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.999185.
# 0 variants removed due to missing genotype data (--geno).
# --hwe: 31 variants removed due to Hardy-Weinberg exact test.
# 7798 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 180616 variants and 3423 people pass filters and QC.
# Note: No phenotypes present.
# Relationship matrix calculation complete.
# --pca: Results saved to final_EUR_ADMIXTURE-PCAS_EUR.eigenval and
# final_EUR_ADMIXTURE-PCAS_EUR.eigenvec .



BFILE="final_AFR_ADMIXTURE"
plink --memory 400000 --threads 24 --bfile ${BFILE} --allow-no-sex --cluster --geno 0.05 --genome --hwe 1e-06 --ld-window-r2 0.2 --maf 0.05 --mds-plot 4 --min 0.2 --nonfounders --pca header --out ${BFILE}-PCAS_AFR
# 188445 variants loaded from .bim file.
# 592 people (0 males, 0 females, 592 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_AFR_ADMIXTURE-PCAS_AFR.nosex .
# Using up to 24 threads (change this with --threads).
# Warning: This run includes BLAS/LAPACK linear algebra operations which
# currently disregard the --threads limit.  If this is problematic, you may want
# to recompile against single-threaded BLAS/LAPACK.
# Before main variant filters, 592 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.999131.
# 23 variants removed due to missing genotype data (--geno).
# --hwe: 56 variants removed due to Hardy-Weinberg exact test.
# 53840 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 134526 variants and 592 people pass filters and QC.
# Note: No phenotypes present.
# Relationship matrix calculation complete.
# --pca: Results saved to final_AFR_ADMIXTURE-PCAS_AFR.eigenval and
# final_AFR_ADMIXTURE-PCAS_AFR.eigenvec .
# IBD calculations complete.
# Finished writing final_AFR_ADMIXTURE-PCAS_AFR.genome .
# Clustering... done.
# Cluster solution written to final_AFR_ADMIXTURE-PCAS_AFR.cluster1 ,
# final_AFR_ADMIXTURE-PCAS_AFR.cluster2 , and
# final_AFR_ADMIXTURE-PCAS_AFR.cluster3 .
# Performing multidimensional scaling analysis (SVD algorithm, 4
# dimensions)... done.
# MDS solution written to final_AFR_ADMIXTURE-PCAS_AFR.mds .



BFILE="final_EAS_ADMIXTURE"
plink --memory 400000 --threads 24 --bfile ${BFILE} --allow-no-sex --cluster --geno 0.05 --genome --hwe 1e-06 --ld-window-r2 0.2 --maf 0.05 --mds-plot 4 --min 0.2 --nonfounders --pca header --out ${BFILE}-PCAS_EAS
# 188445 variants loaded from .bim file.
# 21 people (0 males, 0 females, 21 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to final_EAS_ADMIXTURE-PCAS_EAS.nosex .
# Using up to 24 threads (change this with --threads).
# Warning: This run includes BLAS/LAPACK linear algebra operations which
# currently disregard the --threads limit.  If this is problematic, you may want
# to recompile against single-threaded BLAS/LAPACK.
# Before main variant filters, 21 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.999665.
# 79 variants removed due to missing genotype data (--geno).
# --hwe: 0 variants removed due to Hardy-Weinberg exact test.
# 74651 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 113715 variants and 21 people pass filters and QC.
# Note: No phenotypes present.
# Relationship matrix calculation complete.
# --pca: Results saved to final_EAS_ADMIXTURE-PCAS_EAS.eigenval and
# final_EAS_ADMIXTURE-PCAS_EAS.eigenvec .
# IBD calculations complete.
# Finished writing final_EAS_ADMIXTURE-PCAS_EAS.genome .
# Clustering... done.
# Cluster solution written to final_EAS_ADMIXTURE-PCAS_EAS.cluster1 ,
# final_EAS_ADMIXTURE-PCAS_EAS.cluster2 , and
# final_EAS_ADMIXTURE-PCAS_EAS.cluster3 .
# Performing multidimensional scaling analysis (SVD algorithm, 4
# dimensions)... done.
# MDS solution written to final_EAS_ADMIXTURE-PCAS_EAS.mds .


#####################################################################
## check variants from  SNP eff in Zhaoming and Qin et al's papers ##
#####################################################################
# For this task, I am mostly using R code variants_counts.R
# First create a list of variants to be searched
grep -v '^##' MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr4.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf | cut -d$'\t' -f3| head

rm all_variants.txt
getVars()
{
PATTERN=$1
FILE=$2	
if grep -q ${PATTERN} ${FILE} ; then
echo "found $PATTERN" >> var_search_result.txt
fi
}

export -f getVars

for PATTERN in $(cat zhaoming_variant_sites.txt); do
echo "Doing $PATTERN"
parallel -j22 getVars ${PATTERN} MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr{}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf ::: {1..22}
done

