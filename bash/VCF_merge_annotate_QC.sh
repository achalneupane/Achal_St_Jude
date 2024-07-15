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

## Annotate with the latest version of dbnsfp
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff_v2
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz.tbi .

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=20; \
export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff_v2"; \
export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/12_10_2023/clinvar.vcf.gz"; \
# export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
# export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/dbnsfp/new_dbNSFP4.4a_variant.${CHR}.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
        -P "${CHR}_annotate" \
        -J "${CHR}_ann" \
        -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff_v2/entrypoint_snpEff_annotation_round2.sh"; \
done;


##################################################################
## Helper script to Annotate VCF using snpeff and snpsift tools ##
##################################################################
######################
## Achal Neupane    ##
## Date: 10/10/2023 ##
######################
VERSION="2.0"
#!/usr/bin/bash
# module load gatk/3.7
module load gatk/4.1.8.0
module load vt
module load vcftools
module load bcftools
module load tabix
module load zlib/1.2.5
module load java/13.0.1
module load vep/v108
module load samtools

cd ${WORKDIR}

VCF="${INPUT_VCF}"

MAX_HEADER_LINES=5000
ANNOT_SOURCE="new_${VCF}"; ANNOT_PROJECT="new_${VCF%.*}-annot"

## Adding dbSNP
gatk --java-options "-Xmx16g -XX:ParallelGCThreads=20" VariantAnnotator \
   -R ${REF} \
   -V ${VCF} \
   -L ${VCF} \
   -D /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz \
   -O new_${VCF}

echo "DONE GATK Annotation with dbSNP for ${CHR}" >> annotation_step.txt

## Start annotating
# zcat ${ANNOT_SOURCE} | head -${MAX_HEADER_LINES} | grep "^##" > ${ANNOT_PROJECT}.vcf
zcat ${ANNOT_SOURCE}| grep -v "^##" | cut -f1-8 | awk -F'\t' '!_[$3]++' >> ${ANNOT_PROJECT}.vcf
sed -i 's/\t\*\t/\t<*:DEL>\t/g' ${ANNOT_PROJECT}.vcf
echo "DONE trimming the VCF for ${CHR}" >> annotation_step.txt
## Adding snpEff annotation; these are genome annotations from ENSEMBL, created from GRCh38/hg38 reference genome sequence
${JAVA} ${JAVAOPTS} -jar ${SNPEFF} -v GRCh38.105  ${ANNOT_PROJECT}.vcf > ${ANNOT_PROJECT}-snpeff.vcf
mv snpEff_genes.txt ${CHR}_snpEff_genes.txt; mv snpEff_summary.html ${CHR}_snpEff_summary.html
# rm ${ANNOT_PROJECT}.vcf
echo "DONE SNPeff for ${CHR}" >> annotation_step.txt

# Adding EXaC db
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with EXAC for ${CHR}" >> annotation_step.txt
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf

## Adding clinvar CLNSIG
module load java/13.0.1
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf > ${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.vcf
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with clinvar for ${CHR}" >> annotation_step.txt


## Adding dbNSFP
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -n -v ${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf

# cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_CADD_phred" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_clinvar_clnsig" "dbNSFP_MutationAssessor_pred"      "dbNSFP_MutationTaster_pred"    "dbNSFP_Polyphen2_HDIV_pred"    "dbNSFP_Polyphen2_HVAR_pred"    "dbNSFP_SIFT_pred"      "dbNSFP_LRT_pred"       "CLNSIG" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_raw" "AF_popmax"> ${ANNOTATED%.*}-FIELDS-simple.txt






## Adding latest dbNSFP and gnomad
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff

# gnomAD data
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/gnomAD
for i in {1..22}; do 
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${i}.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${i}.vcf.bgz.tbi 
done


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff_v2
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz.tbi .


for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=2; \
export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff"; \
export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/12_10_2023/clinvar.vcf.gz"; \
# export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
# export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/dbnsfp/new_dbNSFP4.4a_variant.${CHR}.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
 -P "${CHR}_annotate" \
 -q "rhel8_standard" \
 -J "${CHR}_ann" \
 -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
 -n "${THREADS}" \
 -R "rusage[mem=40GB]" \
 "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/latest_dbnsfp.sh"; \
done;


latest_dbnsfp.sh
#!/usr/bin/bash

module load gatk/4.1.8.0
module load vt
module load vcftools
module load bcftools
module load tabix
module load zlib/1.2.5
module load java/13.0.1
module load vep/v108
module load samtools
## Adding latest dbNSFP
cd ${WORKDIR}

# ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -n -v ${WORKDIR}/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf > ${WORKDIR}/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf

ANNOTATED="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf"

## Add gnomad
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/gnomAD/gnomad.genomes.v4.0.sites.${CHR}.vcf.bgz \
${ANNOTATED} > "${WORKDIR}/$(basename ${ANNOTATED} .vcf).gnomAD.vcf"


## Adding clinvar from 12/10/2023
ANNOTATED=MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD.vcf

module load java/13.0.1
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOTATED} > "$(basename ${ANNOTATED} .vcf).gnomAD_clinvar_12_10_2023.vcf"

# ANNOTATED="$(basename ${ANNOTATED} .vcf).gnomAD_clinvar_12_10_2023.vcf"
# cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_CADD_phred" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_clinvar_clnsig" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationTaster_pred" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_SIFT_pred" "dbNSFP_LRT_pred" "CLNSIG" "dbNSFP_gnomAD_genomes_AF" "dbNSFP_gnomAD_genomes_NFE_AF" "dbNSFP_gnomAD_genomes_AMI_AF" "dbNSFP_gnomAD_genomes_ASJ_AF" "dbNSFP_gnomAD_genomes_MID_AF" "dbNSFP_gnomAD_genomes_FIN_AF" "dbNSFP_gnomAD_genomes_SAS_AF" "dbNSFP_gnomAD_genomes_EAS_AF" "dbNSFP_gnomAD_genomes_AFR_AF" "dbNSFP_gnomAD_genomes_POPMAX_AF" "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_exomes_NFE_AF" "dbNSFP_gnomAD_exomes_ASJ_AF" "dbNSFP_gnomAD_exomes_FIN_AF" "dbNSFP_gnomAD_exomes_SAS_AF" "dbNSFP_gnomAD_exomes_EAS_AF" "dbNSFP_gnomAD_exomes_AFR_AF" "dbNSFP_gnomAD_exomes_POPMAX_AF" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_fin" "AF_raw" "AF" "AF_grpmax" "AF_joint_nfe" "AF_joint_afr" "AF_joint_eas" "AF_joint_sas" "AF_joint_fin" "AF_joint_raw" "AF_joint" "AF_grpmax_joint" > ${WORKDIR}/${ANNOTATED%.*}-clinvar_12_10_2023_FIELDS-simple2.txt
ANNOTATED="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD.gnomAD_clinvar_12_10_2023.vcf"
OUTFILE="${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt"
cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "dbNSFP_gnomAD_genomes_AF" "dbNSFP_gnomAD_genomes_NFE_AF" "dbNSFP_gnomAD_genomes_AMI_AF" "dbNSFP_gnomAD_genomes_ASJ_AF" "dbNSFP_gnomAD_genomes_MID_AF" "dbNSFP_gnomAD_genomes_FIN_AF" "dbNSFP_gnomAD_genomes_SAS_AF" "dbNSFP_gnomAD_genomes_EAS_AF" "dbNSFP_gnomAD_genomes_AFR_AF" "dbNSFP_gnomAD_genomes_POPMAX_AF" "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_exomes_NFE_AF" "dbNSFP_gnomAD_exomes_ASJ_AF" "dbNSFP_gnomAD_exomes_FIN_AF" "dbNSFP_gnomAD_exomes_SAS_AF" "dbNSFP_gnomAD_exomes_EAS_AF" "dbNSFP_gnomAD_exomes_AFR_AF" "dbNSFP_gnomAD_exomes_POPMAX_AF" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_fin" "AF_raw" "AF" "AF_grpmax" "AF_joint_nfe" "AF_joint_afr" "AF_joint_eas" "AF_joint_sas" "AF_joint_fin" "AF_joint_raw" "AF_joint" "AF_grpmax_joint" > ${WORKDIR}/${OUTFILE}



awk -F'\t' '$48 < 0.01 && $54 < 0.01 && $48 != "." && $54 != "."' test.field.txt | less -S

################
## 07/12/2023 ##
################
## Adding annotations with loftee latest dbnsfp
mkdir snpEff_v2

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff_v2
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz.tbi .


for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=2; \
export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff"; \
export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/12_10_2023/clinvar.vcf.gz"; \
# export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
# export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/dbnsfp/new_dbNSFP4.4a_variant.${CHR}.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
 -P "${CHR}_annotate" \
 -q "rhel8_standard" \
 -J "${CHR}_ann" \
 -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
 -n "${THREADS}" \
 -R "rusage[mem=40GB]" \
 "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/latest_dbnsfp.sh"; \
done;


latest_dbnsfp.sh
#!/usr/bin/bash

module load gatk/4.1.8.0
module load vt
module load vcftools
module load bcftools
module load tabix
module load zlib/1.2.5
module load java/13.0.1
module load vep/v108
module load samtools
## Adding latest dbNSFP
cd ${WORKDIR}

# ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -n -v ${WORKDIR}/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf > ${WORKDIR}/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf

ANNOTATED="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf"

## Add gnomad
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/gnomAD/gnomad.genomes.v4.0.sites.${CHR}.vcf.bgz \
${ANNOTATED} > "${WORKDIR}/$(basename ${ANNOTATED} .vcf).gnomAD.vcf"


## Adding clinvar from 12/10/2023
ANNOTATED=MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD.vcf

module load java/13.0.1
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOTATED} > "$(basename ${ANNOTATED} .vcf).gnomAD_clinvar_12_10_2023.vcf"

# ANNOTATED="$(basename ${ANNOTATED} .vcf).gnomAD_clinvar_12_10_2023.vcf"
# cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_CADD_phred" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_clinvar_clnsig" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationTaster_pred" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_SIFT_pred" "dbNSFP_LRT_pred" "CLNSIG" "dbNSFP_gnomAD_genomes_AF" "dbNSFP_gnomAD_genomes_NFE_AF" "dbNSFP_gnomAD_genomes_AMI_AF" "dbNSFP_gnomAD_genomes_ASJ_AF" "dbNSFP_gnomAD_genomes_MID_AF" "dbNSFP_gnomAD_genomes_FIN_AF" "dbNSFP_gnomAD_genomes_SAS_AF" "dbNSFP_gnomAD_genomes_EAS_AF" "dbNSFP_gnomAD_genomes_AFR_AF" "dbNSFP_gnomAD_genomes_POPMAX_AF" "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_exomes_NFE_AF" "dbNSFP_gnomAD_exomes_ASJ_AF" "dbNSFP_gnomAD_exomes_FIN_AF" "dbNSFP_gnomAD_exomes_SAS_AF" "dbNSFP_gnomAD_exomes_EAS_AF" "dbNSFP_gnomAD_exomes_AFR_AF" "dbNSFP_gnomAD_exomes_POPMAX_AF" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_fin" "AF_raw" "AF" "AF_grpmax" "AF_joint_nfe" "AF_joint_afr" "AF_joint_eas" "AF_joint_sas" "AF_joint_fin" "AF_joint_raw" "AF_joint" "AF_grpmax_joint" > ${WORKDIR}/${ANNOTATED%.*}-clinvar_12_10_2023_FIELDS-simple2.txt

# ANNOTATED="$(basename ${ANNOTATED} .vcf).gnomAD_clinvar_12_10_2023.vcf"
ANNOTATED="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD.gnomAD_clinvar_12_10_2023.vcf"
OUTFILE="${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple3.txt"
cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "dbNSFP_gnomAD_genomes_AF" "dbNSFP_gnomAD_genomes_NFE_AF" "dbNSFP_gnomAD_genomes_AMI_AF" "dbNSFP_gnomAD_genomes_ASJ_AF" "dbNSFP_gnomAD_genomes_MID_AF" "dbNSFP_gnomAD_genomes_FIN_AF" "dbNSFP_gnomAD_genomes_SAS_AF" "dbNSFP_gnomAD_genomes_EAS_AF" "dbNSFP_gnomAD_genomes_AFR_AF" "dbNSFP_gnomAD_genomes_POPMAX_AF" "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_exomes_NFE_AF" "dbNSFP_gnomAD_exomes_ASJ_AF" "dbNSFP_gnomAD_exomes_FIN_AF" "dbNSFP_gnomAD_exomes_SAS_AF" "dbNSFP_gnomAD_exomes_EAS_AF" "dbNSFP_gnomAD_exomes_AFR_AF" "dbNSFP_gnomAD_exomes_POPMAX_AF" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_fin" "AF_raw" "AF" "AF_grpmax" "AF_joint_nfe" "AF_joint_afr" "AF_joint_eas" "AF_joint_sas" "AF_joint_fin" "AF_joint_raw" "AF_joint" "AF_grpmax_joint" > ${WORKDIR}/${$OUTFILE}



## I was not able to split the INFO field for Loftee, so keeping them in a separate folder
mkdir loftee
## Loftee
vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
--offline --everything --force_overwrite \
--assembly GRCh38 \
--dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ \
--fasta /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
-i ${ANNOTATED} \
--tab -o ${WORKDIR}/loftee/${CHR}.preQC_biallelic_renamed_ID.loftee.tsv \
--plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:false


echo "DONE for ${CHR}" >> annotation_step.txt



#########

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/loftee
ls chr*.preQC_biallelic_renamed_ID.loftee.tsv| sort -V| xargs cat| grep -v ^##| awk -F'\t' 'NR==1 || ($81 ~ /HC/)' > Loftee_HC_all_chr.txt
awk 'NR > 1 {print $1}' Loftee_HC_all_chr.txt| sort -V | uniq > loftee_unique_HC_variants.txt
## Extract annotated file lines that match HC
for i in {1..22}; do 
export CHR="chr${i}"; 
echo "Working on $CHR"; 
head -1 ../${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > ${CHR}_loftee_HC_variant_extracted_v4.txt
grep -F -f loftee_unique_HC_variants.txt ../${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt >> ${CHR}_loftee_HC_variant_extracted_v4.txt
done

## Extract missense variants
for i in {1..22}; do 
export CHR="chr${i}"; 
echo "Working on $CHR"; 
head -1 ${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > ${CHR}_missense_variant_v4.txt
grep missense ${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt >> ${CHR}_missense_variant_v4.txt
done



## P/LP clinvar
cut -d$'\t' -f210 chr*.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt| sort -V | uniq
awk -F'\t' 'NR==1 || ($210 ~ /Pathogenic|Likely_pathogenic/) && $210 !~ /Conflicting|Benign/' \
chr*.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > all_new_clinvar_P_LP.txt

# (base) [aneupane@splprhpc07 snpEff_round3]$ wc -l all_new_clinvar_P_LP.txt
# 99787 all_new_clinvar_P_LP.txt

# extract unique
awk -F'\t' '!seen[$3]++' all_new_clinvar_P_LP.txt > all_new_clinvar_P_LP_unique.txt







#####################################
## Annovar and intervar Annotation ##
#####################################
# First check avblist hg38_avdblist.txt for the latest releases
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar
# perl annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg38 /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/annovar/humandb/

# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar 1000g2015aug humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03nontcga humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar esp6500siv2_all humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_exome humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar revel humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar intervar_20180118 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbscsnv11 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic70 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar nci60 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20220320 humandb/ 


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
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/12_10_2023/clinvar.vcf.gz"; \
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


####################################################################
## check variants from  SNPeff in Zhaoming and Qin et al's papers ##
####################################################################
# For this task, I am mostly using R code variants_counts.R
# First create a list of variants to be searched

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/variants_in_VCF
ln -s ../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf .

module load parallel
getVars()
{
VCF="$1"
CHR="$2"
grep -v "^##" ${VCF} | cut -f3 | cut -d";" -f1 > CHR${CHR}_Vars.vcf
}

export -f getVars
parallel -j22 getVars MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr{}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf {} ::: {1..22}


# Then run R script variants_counts.R

##################
## LoF variants ##
##################
# Qin et al: frameshift      missense      nonsense    proteinDel        splice splice_region 
# 196            70           175             1            74             4 

## Annotation of interest
# frameshift_variant, start_lost, stop_gained, splice, ^splice_region_variant$
head -1 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt > LoF_variants.txt
cat $(ls *.dbSNP155-FIELDS-simple.txt| sort -V)| egrep 'frameshift_variant|start_lost|stop_gained|splice|splice_region_variant' >> LoF_variants.txt
awk '{ print $3 }' LoF_variants.txt| cut -d";" -f1 > LoF_variants_ID.txt

## Now I am checking the concordance for Indels from the Qin et al and Zhaoming et al in our VCF file
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/attr_fraction
zhaoming_et_al_variants_bed.txt
awk '{if($2 != $3) print}' zhaoming_et_al_variants_bed.txt > zhaoming_et_al_variants_INDEL.bed
awk '{if($2 != $3) print}' qin_et_al_variants_bed.txt > qin_et_al_variants_INDEL.bed


ln -s ../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf .

tabIndexed()
{
module load bcftools
module load tabix	
CHR="$1"
VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf"
bgzip -c ${VCF} > ${VCF}.gz
tabix -p vcf ${VCF}.gz
}

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/annotated_indexed_vcf
export -f tabIndexed
parallel -j22 tabIndexed {} ::: {1..22}


## Extract all INDELs from the annotated VCF; These VCF needs to be tabix indexed
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/annotated_indexed_vcf
ln -s ../zhaoming_et_al_variants_INDEL.bed .
ln -s ../qin_et_al_variants_INDEL.bed .

#!/bin/bash

# Extract genotypes from the original GATK VCF (provided by Comp. Bio. department) before any QC/hard filtering
module load bcftools

############################################################################################################################################
## Now, let's search all the INDELS in Zhaoming and Qin et al in our annotated VCF files. We are checking indels +-1bp up and down stream ##
############################################################################################################################################
################################
## INDELS from Zhaoming et al ##
################################
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/annotated_indexed_vcf

# ## Tabix
# # line="chr1:17028712:_:C       chr1    17028711        17028714        chr1:17028711-17028714"
# zgrep "^#CHROM" MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr22.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz|head -1 | awk '{ print "KEY.varID\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> zhaoming_et_al_variants_INDEL.bed.out
# LIST=zhaoming_et_al_variants_INDEL.bed
# # tabixSearch()
# # {
# # #!/bin/bash
# # module load tabix
# # CHR="$1"
# # LIST="$2"
# for CHR in {1..22}; do
# VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz"
# for line in $(cat ${LIST}); do
# FOUND_LINE=()	
# query="$(echo ${line}| awk '{ print $5 }')"
# SNPId="$(echo ${line}| awk '{ print $1 }')"
# IFS=$'\n'
# FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
# if [ -n "${FOUND_LINE}" ]; then
# for each in "${FOUND_LINE[@]}"
# do
# echo -e "${SNPId}\t${each}" >> ${LIST}.out
# done
# fi
# done
# done
# # }



# ## Search in +/- 5 bases
# LIST=zhaoming_et_al_variants_INDEL.bed
# zgrep "^#CHROM" MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr22.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz|head -1 | awk '{ print "KEY.varID\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> ${LIST}_plus_minu10bps.out
# for CHR in {1..22}; do
# VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz"
# for line in $(cat ${LIST}); do
# FOUND_LINE=()	
# query="$(echo ${line}| awk '{ print $6 }')"
# SNPId="$(echo ${line}| awk '{ print $1 }')"
# IFS=$'\n'
# FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
# if [ -n "${FOUND_LINE}" ]; then
# for each in "${FOUND_LINE[@]}"
# do
# echo -e "${SNPId}\t${each}" >> ${LIST}_plus_minu10bps.out
# done
# fi
# done
# done



# # export -f tabixSearch
# # export LIST="zhaoming_et_al_variants_INDEL.bed"
# # parallel -j22 tabixSearch {} ${LIST} ::: {1..22}

# # ## Bcftools
# # zgrep "^#CHROM" MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr22.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz|head -1| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' > sjlife_zhaoming.vcf
# # for chr in {1..22}; do
# # bcftools view -Ov MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz -R zhaoming_et_al_variants_INDEL.bed | grep -v '^#'| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' >> sjlife_zhaoming.vcf
# # done


# ###########################
# ## INDELS from Qin et al ##
# ###########################
# ## Tabix
# zgrep "^#CHROM" MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr22.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz|head -1 | awk '{ print "KEY.varID\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> qin_et_al_variants_INDEL.bed.out
# LIST=qin_et_al_variants_INDEL.bed
# for CHR in {1..22}; do
# VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz"
# for line in $(cat ${LIST}); do
# FOUND_LINE=()	
# query="$(echo ${line}| awk '{ print $5 }')"
# SNPId="$(echo ${line}| awk '{ print $1 }')"
# IFS=$'\n'
# FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
# if [ -n "${FOUND_LINE}" ]; then
# for each in "${FOUND_LINE[@]}"
# do
# echo -e "${SNPId}\t${each}" >> ${LIST}.out
# done
# fi
# done
# done


# ## Search in +/- 5 bases
# LIST=qin_et_al_variants_INDEL.bed
# zgrep "^#CHROM" MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr22.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz|head -1 | awk '{ print "KEY.varID\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> ${LIST}_plus_minu10bps.out
# for CHR in {1..22}; do
# VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz"
# for line in $(cat ${LIST}); do
# FOUND_LINE=()	
# query="$(echo ${line}| awk '{ print $6 }')"
# SNPId="$(echo ${line}| awk '{ print $1 }')"
# IFS=$'\n'
# FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
# if [ -n "${FOUND_LINE}" ]; then
# for each in "${FOUND_LINE[@]}"
# do
# echo -e "${SNPId}\t${each}" >> ${LIST}_plus_minu10bps.out
# done
# fi
# done
# done

## saved this bit of code as search_indel.sh

# # Bcftools
# zgrep "^#CHROM" MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr22.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz|head -1| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' > sjlife_qin.vcf
# for chr in {1..22}; do
# bcftools view -Ov MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz -R qin_et_al_variants_INDEL.bed | grep -v '^#'| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' >> sjlife_qin.vcf
# done

## Now, I will manually check each indel for exact match.
# From Zhaoming's list of INDELs, I was able to find only 70/133 INDELs in our VCF
# From Qin's list of INDELs, I was able to find only 82/181 INDELs in our VCF
# This is saved as zhaoming_and_qin_et_al_variant_INDEL_comparison_in_SJLIFE.xlxs on /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/annotated_indexed_vcf


# Also checking these in Yadav's VCF files
VCF="sjlife_1_zhaoming.vcf"
zcat $VCF | bgzip -c > ${VCF}.gz
tabix -p vcf ${VCF}.gz


for f in *_v2.vcf; do
	echo $f
	cp $f ${f}.gz
done

for f in *_v2.vcf.gz; do
	echo $f
bcftools sort -Oz $f -o sorted_${f}
done

for f in sorted_*v2.vcf.gz; do
	echo $f
bcftools index -f -t --threads 4 ${f}
done


## Edit Zhaoming vcf
bcftools merge --threads 4 sorted_sjlife_1_zhaoming_v2.vcf.gz sorted_sjlife_2_zhaoming_v2.vcf.gz -0 -Oz -o MERGED_sorted_sjlife_1_2_zhaoming_v2.vcf.gz
bcftools index -f -t --threads 4 MERGED_sorted_sjlife_1_2_zhaoming_v2.vcf.gz
## Convert to biallelic
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa MERGED_sorted_sjlife_1_2_zhaoming_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_zhaoming_v2.vcf.gz
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' MERGED_biallelic_sorted_sjlife_1_2_zhaoming_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz
bcftools index -f -t --threads 4 MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz


## Edit qin vcf
bcftools merge --threads 4 sorted_sjlife_1_qin_v2.vcf.gz sorted_sjlife_2_qin_v2.vcf.gz -0 -Oz -o MERGED_sorted_sjlife_1_2_qin_v2.vcf.gz
bcftools index -f -t --threads 4 MERGED_sorted_sjlife_1_2_qin_v2.vcf.gz
## Convert to biallelic
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa MERGED_sorted_sjlife_1_2_qin_v2.vcf.gz Oz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_qin_v2.vcf.gz
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' MERGED_biallelic_sorted_sjlife_1_2_qin_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_qin_ID_edited.vcf.gz
bcftools index -f -t --threads 4 MERGED_biallelic_sorted_sjlife_1_2_qin_ID_edited.vcf.gz



## Search in +/- 10 bases in Yadav's VCF (VCF prior to QC)
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/attr_fraction
# ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/annotated_indexed_vcf/zhaoming_et_al_variants_INDEL.bed .
# ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/annotated_indexed_vcf/qin_et_al_variants_INDEL.bed .

# LIST=zhaoming_et_al_variants_INDEL.bed
# VCF="MERGED_sorted_sjlife_1_2_zhaoming.vcf.gz"
# zgrep "^#CHROM" ${VCF}|head -1 | awk '{ print "KEY.varID\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> ${LIST}_plus_minu10bps.out
# for line in $(cat ${LIST}); do
# FOUND_LINE=()	
# query="$(echo ${line}| awk '{ print $6 }')"
# SNPId="$(echo ${line}| awk '{ print $1 }')"
# IFS=$'\n'
# FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
# if [ -n "${FOUND_LINE}" ]; then
# for each in "${FOUND_LINE[@]}"
# do
# echo -e "${SNPId}\t${each}" >> ${LIST}_plus_minu10bps.out
# done
# else
# echo -e "${SNPId}\tNA\tNA\tNA\tNA" >> ${LIST}_plus_minu10bps.out
# fi
# done


# LIST=qin_et_al_variants_INDEL.bed
# VCF="MERGED_sorted_sjlife_1_2_qin.vcf.gz"
# zgrep "^#CHROM" ${VCF}|head -1 | awk '{ print "KEY.varID\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> ${LIST}_plus_minu10bps.out
# for line in $(cat ${LIST}); do
# FOUND_LINE=()	
# query="$(echo ${line}| awk '{ print $6 }')"
# SNPId="$(echo ${line}| awk '{ print $1 }')"
# IFS=$'\n'
# FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
# if [ -n "${FOUND_LINE}" ]; then
# for each in "${FOUND_LINE[@]}"
# do
# echo -e "${SNPId}\t${each}" >> ${LIST}_plus_minu10bps.out
# done
# fi
# done


# Search in +/- 10 bases in Yadav's VCF (VCF prior to QC)
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/attr_fraction
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/zhaoming_et_al_SNV_INDEL_Search_list_V2.bed .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/qin_et_al_SNV_INDEL_Search_list_V2.bed .

LIST="zhaoming_et_al_SNV_INDEL_Search_list_V2.bed"
VCF="MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz"
zgrep "^#CHROM" ${VCF}|head -1 | awk '{ print "KEY.varID\tVAR_TYPE\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> ${LIST}_V2_results.out
for line in $(cat ${LIST}); do
FOUND_LINE=()	
query="$(echo ${line}| awk '{ print $6 }')"
VAR_TYPE="$(echo ${line}| awk '{ print $7 }')"
SNPId="$(echo ${line}| awk '{ print $1 }')"
IFS=$'\n'
FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
if [ -n "${FOUND_LINE}" ]; then
for each in "${FOUND_LINE[@]}"
do
echo -e "${SNPId}\t${VAR_TYPE}\t${each}" >> ${LIST}_V2_results.out
done
else
echo -e "${SNPId}\t${VAR_TYPE}\tNA\tNA\tNA\tNA\tNA" >> ${LIST}_V2_results.out
fi
done
sed -i "s/\r//g"  ${LIST}_V2_results.out


LIST="qin_et_al_SNV_INDEL_Search_list_V2.bed"
VCF="MERGED_biallelic_sorted_sjlife_1_2_qin_ID_edited.vcf.gz"
zgrep "^#CHROM" ${VCF}|head -1 | awk '{ print "KEY.varID\tVAR_TYPE\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> ${LIST}_V2_results.out
for line in $(cat ${LIST}); do
FOUND_LINE=()	
query="$(echo ${line}| awk '{ print $6 }')"
VAR_TYPE="$(echo ${line}| awk '{ print $7 }')"
SNPId="$(echo ${line}| awk '{ print $1 }')"
IFS=$'\n'
FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
if [ -n "${FOUND_LINE}" ]; then
for each in "${FOUND_LINE[@]}"
do
echo -e "${SNPId}\t${VAR_TYPE}\t${each}" >> ${LIST}_V2_results.out
done
else
echo -e "${SNPId}\t${VAR_TYPE}\tNA\tNA\tNA\tNA\tNA" >> ${LIST}_V2_results.out
fi
done
sed -i "s/\r//g"  ${LIST}_V2_results.out

@@@@@ Open these two output files on Excel spreadsheet `Zhaoming_and_qin_in_VCF_prior_to_any_QC.xlxs` on Z:\ResearchHome\Groups\sapkogrp\projects\Genomics\common\attr_fraction
@@@@@ I have checked the concordance of these variants manually.




#####################################################
## Now search these variants in the VCF I prepared ##
#####################################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/annotated_indexed_vcf

# Let's merge the vcf.gz 
module load bcftools
FILES=$(ls MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz| sort -V)
bcftools concat ${FILES} -Oz -o MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr.ALL.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz 
bcftools index -f -t --threads 20 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr.ALL.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz 


ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/zhaoming_et_al_SNV_INDEL_Search_list_V2.bed .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/qin_et_al_SNV_INDEL_Search_list_V2.bed .

LIST="zhaoming_et_al_SNV_INDEL_Search_list_V2.bed"
VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr.ALL.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz"
zcat ${VCF} | head -10000 | zgrep "^#CHROM"|head -1 | awk '{ print "KEY.varID\tVAR_TYPE\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> POST_QC_VCF_${LIST}_V2_results.out
for line in $(cat ${LIST}); do
FOUND_LINE=()	
query="$(echo ${line}| awk '{ print $6 }')"
VAR_TYPE="$(echo ${line}| awk '{ print $7 }')"
SNPId="$(echo ${line}| awk '{ print $1 }')"
IFS=$'\n'
FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
if [ -n "${FOUND_LINE}" ]; then
for each in "${FOUND_LINE[@]}"
do
echo -e "${SNPId}\t${VAR_TYPE}\t${each}" >> POST_QC_VCF_${LIST}_V2_results.out
done
else
echo -e "${SNPId}\t${VAR_TYPE}\tNA\tNA\tNA\tNA\tNA" >> POST_QC_VCF_${LIST}_V2_results.out
fi
done
sed -i "s/\r//g"  POST_QC_VCF_${LIST}_V2_results.out


LIST="qin_et_al_SNV_INDEL_Search_list_V2.bed"
VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr.ALL.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.gz"
zcat ${VCF} | head -10000 | zgrep "^#CHROM"|head -1 | awk '{ print "KEY.varID\tVAR_TYPE\t"$1"\t"$2"\t"$3"\t"$4"\t"$5 }'> POST_QC_VCF_${LIST}_V2_results.out
for line in $(cat ${LIST}); do
FOUND_LINE=()	
query="$(echo ${line}| awk '{ print $6 }')"
VAR_TYPE="$(echo ${line}| awk '{ print $7 }')"
SNPId="$(echo ${line}| awk '{ print $1 }')"
IFS=$'\n'
FOUND_LINE=( $(tabix ${VCF} ${query}| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' ) )
if [ -n "${FOUND_LINE}" ]; then
for each in "${FOUND_LINE[@]}"
do
echo -e "${SNPId}\t${VAR_TYPE}\t${each}" >> POST_QC_VCF_${LIST}_V2_results.out
done
else
echo -e "${SNPId}\t${VAR_TYPE}\tNA\tNA\tNA\tNA\tNA" >> POST_QC_VCF_${LIST}_V2_results.out
fi
done
sed -i "s/\r//g"  POST_QC_VCF_${LIST}_V2_results.out


# Now use R code: variant_counts.R

# #####################
# ## Filter variants ##
# #####################
# CLINVAR <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "Clinvar",]
# CLINVAR <- CLINVAR[!duplicated(CLINVAR$KEY.varID),]

# MetaSVM <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "MetaSVM",]
# MetaSVM <- MetaSVM[!duplicated(MetaSVM$KEY.varID),]

# LOF.VCF <- fread("LoF_variants.txt", header = T, sep = "\t")
# LOF.VCF$PRED_TYPE <- "LoF"
# LOF.VCF$KEY.pos <- paste(LOF.VCF$CHROM, LOF.VCF$POS, sep = ":")        
# LOF.VCF$KEY.varID <- paste(LOF.VCF$CHROM, LOF.VCF$POS, LOF.VCF$REF, LOF.VCF$ALT, sep = ":")     

# non.intron.LoF <- LOF.VCF[grepl("splice_region_variant&intron_variant", LOF.VCF$`ANN[*].EFFECT`),]
# sum(unique(zhaoming.SNV$KEY.varID) %in% non.intron.LoF$KEY.varID)
# tt <- non.intron.LoF[non.intron.LoF$KEY.varID %in% unique(zhaoming.SNV$KEY.varID),]

# # LOF.VCF <- LOF.VCF[!duplicated(LOF.VCF$KEY.varID),]

# #####################################################################
# #####################################################################
# #####################################################################
# #####################################################################

# ##################################################
# ## Now, check these in original VCFs from Yadav ##
# ##################################################
# ## Zhaoming
# zhaoming_in_vcf <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/Zhaoming_in_VCF_prior_to_any_QC.txt", sep = "\t", header = T, check.names = F)
# zhaoming_in_vcf <- zhaoming_in_vcf[1:10]
# head(zhaoming_in_vcf)
# # zhaoming_in_vcf <- zhaoming_in_vcf[zhaoming_in_vcf$Match_per_var == "Y",]

# zhaoming_in_vcf$MATCHED.IN.CLINVAR <- ""
# zhaoming_in_vcf$MATCHED.IN.MetaSVM <- ""
# zhaoming_in_vcf$MATCHED.IN.LOF <- ""
# zhaoming_in_vcf$MATCHED.IN.ANY.ANNOTATION <- ""
# zhaoming_in_vcf$MATCHED.IN.CLINVAR[zhaoming_in_vcf$ID %in% CLINVAR$KEY.varID] <- "Y"
# zhaoming_in_vcf$MATCHED.IN.MetaSVM[zhaoming_in_vcf$ID %in% MetaSVM$KEY.varID] <- "Y"
# zhaoming_in_vcf$MATCHED.IN.LOF[zhaoming_in_vcf$ID %in% LOF.VCF$KEY.varID] <- "Y"
# zhaoming_in_vcf$MATCHED.IN.ANY.ANNOTATION[zhaoming_in_vcf$ID %in% unique(c(CLINVAR$KEY.varID, MetaSVM$KEY.varID, LOF.VCF$KEY.varID))] <- "Y"

# # Now if zhaoming_in_vcf$Match_YN is not Y, then make all columns Not Y
# zhaoming_in_vcf[!zhaoming_in_vcf$Match_YN == "Y",c("MATCHED.IN.CLINVAR", "MATCHED.IN.MetaSVM", "MATCHED.IN.LOF", "MATCHED.IN.ANY.ANNOTATION")] <- ""
# zhaoming_in_vcf$notes <- ""
# ## Note: Match_per_var is where the variants in VCF match exactly (YES/NO). It is a unique Y/N value per variant whereas Match_YN could have duplicate YES/or NOs
# library(dplyr)
# library(tidyr)
# # Return column names in a column where the value matches a given string
# zhaoming_in_vcf <- zhaoming_in_vcf %>% 
#   mutate(across(c(MATCHED.IN.CLINVAR, MATCHED.IN.MetaSVM, MATCHED.IN.LOF), ~case_when(. == "Y" ~ cur_column()), .names = 'new_{col}')) %>%
#   unite(notes, starts_with('new'), na.rm = TRUE, sep = ';')
# zhaoming_in_vcf$notes <- gsub("MATCHED.IN.", "", zhaoming_in_vcf$notes)

# zhaoming_in_vcf$notes [zhaoming_in_vcf$Match_per_var =="Y" & !zhaoming_in_vcf$MATCHED.IN.ANY.ANNOTATION == "Y" ] <- "QC_DROPPED"
# zhaoming_in_vcf$notes [is.na(zhaoming_in_vcf$ID)] <- "NOT_CALLED_BY_GATK"

# write.table(zhaoming_in_vcf, " Zhaoming_and_qin_in_VCF_prior_to_any_QC_Final_list.txt", row.names = F, col.names = F, quote = F, sep = "\t")

#################################################################
## Checking if any of the variants not called by GATK present: ##
#################################################################
# ## Zhaoming:
# for var in chr12:25245328:C:G chr13:48473359:G:A chr17:7674230:C:T chr12:111418878:G:A

# arr_variable=(chr12:25245328-25245328 chr13:48473359-48473359 chr17:7674230-7674230 chr12:111418878-111418878)

# ## now loop through the above array
# VCF="sorted_sjlife_1_zhaoming_v2.vcf.gz"
# VCF="sorted_sjlife_2_zhaoming_v2.vcf.gz"

# for i in "${arr_variable[@]}"
# do
# tabix ${VCF} $i| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }'
# done

# arr_variable=(chr5:75584794-75584794)
# VCF="sorted_sjlife_1_qin_v2.vcf.gz"
# VCF="sorted_sjlife_2_qin_v2.vcf.gz"

# for i in "${arr_variable[@]}"
# do
# tabix ${VCF} $i| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }'
# done

# ## None were found in individual VCFs

module load vcftools
TYPE="zhaoming"
# TYPE="qin"
VCF="MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.vcf.gz"


vcftools --gzvcf ${VCF} \
--minGQ 5 \
--min-meanDP 5 \
--recode \
--stdout \
| bgzip > MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq5.no.filter.pass.vcf.gz


vcftools --gzvcf ${VCF} \
--keep-filtered PASS \
--minGQ 5 \
--min-meanDP 5 \
--recode \
--stdout \
| bgzip > MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq5.pass.vcf.gz


vcftools --gzvcf ${VCF} \
--keep-filtered PASS \
--minGQ 10 \
--min-meanDP 5 \
--recode \
--stdout \
| bgzip > MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq10.pass.vcf.gz

vcftools --gzvcf ${VCF} \
--keep-filtered PASS \
--minGQ 15 \
--min-meanDP 10 \
--recode \
--stdout \
| bgzip > MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq15.pass.vcf.gz

vcftools --gzvcf ${VCF} \
--keep-filtered PASS \
--minGQ 20 \
--min-meanDP 10 \
--recode \
--stdout \
| bgzip > MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq20.pass.vcf.gz


vcftools --gzvcf ${VCF} \
--minGQ 20 \
--min-meanDP 10 \
--recode \
--stdout \
| bgzip > MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq20.no.filterpass.vcf.gz


# Now let's check the vars determined to be dropped in qin VCF are really dropped after applying this QC
# echo $(cat vars_dropped_in_QC_qin.txt| paste -s -d"|")
# wc -l vars_dropped_in_QC_qin.txt
# 106 vars
TYPE="qin"
var="chr14:75047852:T:TAACTCGCCCATAACTA|chr13:32339287:AAAAG:A|chr16:1775960:TTCCCC:T|chr2:127272934:CT:C|chr6:30708301:ACT:A|chr8:42355545:CATCATCAGCGAATTGGGCTGAAGTAAG:C|chr6:31759819:GCCAGGTGCTGGCA:G|chr15:30905539:GGTTGAAAAACGTGAGGCAT:G|chr17:43070931:TTCTTCTGGGGTCAGGCCAG:T|chr13:102872217:C:CA|chr3:121522031:TAC:T|chr1:27894301:CCA:C|chr17:58732547:CAGT:C|chr9:35074172:GGACGGATCCA:G|chr13:32362626:GCCTTT:G|chr7:7639125:AT:A|chr8:89984550:CA:C|chr19:49862556:GAA:G|chr15:89333601:G:GCTA|chr3:14147980:GCAGACGATGTATC:G|chr12:21475808:AAAAC:A|chr11:65860844:GAGGCGACCCGCAGC:G|chr13:102862141:CAG:C|chr12:21471424:ATACT:A|chr6:30707797:ACT:A|chr6:30708067:G:GT|chr15:89324187:CCT:C|chr13:32337160:TAAAC:T|chr19:45364832:CCTCA:C|chr2:68490198:CAA:C|chr14:45159189:C:CA|chr14:75030646:CT:C|chr7:44073203:A:AG|chr17:43094316:TGATTCAGACTCCCCATCATGTGAGTCATCAGAACCTAACA:T|chr10:49471091:CTA:C|chr10:101579579:AG:A|chr17:35101325:AGCCTCCC:A|chr2:127286924:C:CACTGCT|chr4:83462857:CAAAG:C|chr6:35457937:G:GC|chr15:30905892:T:TA|chr5:132609339:CAAAG:C|chr13:32379464:AC:A|chr17:43124027:ACT:A|chr10:49461381:CCT:C|chr1:241861444:GA:G|chr15:90765412:GATGTAAGTT:G|chr14:45189222:AC:A|chr2:127289348:CAA:C|chr8:11786115:C:CAGCA|chr7:44080814:T:TG|chr5:80672319:C:CTG|chr17:43092848:GTT:G|chr5:132604019:CAA:C|chr16:3589547:CCT:C|chr10:49482711:T:TC|chr1:75815039:A:AT|chr11:65862211:AAT:A|chr10:49472935:CTGGGTCATAGA:C|chr17:50376072:T:TA|chr9:35079169:G:GC|chr1:46273464:AG:A|chr1:46260006:TC:T|chr13:32337160:TA:T|chr1:2201225:CCT:C|chr4:177341536:GGT:G|chr10:14936586:TCTTC:T|chr14:20352261:A:AC|chr1:241878811:ATG:A|chr13:32339482:TTATG:T|chr3:129433215:GA:G|chr13:32339062:CAG:C|chr7:5997348:CG:C|chr9:97689554:CA:C|chr14:45137187:GTGTT:G|chr19:48137531:AG:A|chr16:13926709:AAG:A|chr16:89740841:GAGA:G|chr11:22625192:CCGCTATCACCTTCAGGAAGTTGTTCTGAGGCAAG:C|chr17:34983545:C:CAT|chr16:23635228:AC:A|chr5:176908638:TG:T|chr5:132589652:ACT:A|chr1:226365988:G:GT|chr10:101579754:AC:A|chr19:43560989:TC:T|chr17:61715948:T:TA|chr17:34983189:AG:A|chr3:14156338:GAGTAGACCGC:G|chr3:121489730:TA:T|chr5:132557417:T:TA|chr8:94387020:AG:A|chr2:47799609:CAAAG:C|chr5:75569452:T:TG|chr15:43457065:C:CT|chr15:75349042:G:GC|chr1:46274204:CA:C|chr14:60734947:ACT:A|chr2:68526220:AAG:A|chr19:49862415:GA:G|chr7:44074235:T:G|chr11:47234610:C:T|chr11:65863848:C:T|chr11:65861036:C:T|chr11:61331693:T:G"
echo -e "${TYPE} vars from VCF with no filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.vcf.gz| wc -l)"  > variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp5, gq5 NOpass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq5.no.filter.pass.vcf.gz| wc -l)"  >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp5, gq5 pass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq5.pass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp5, gq10 pass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq10.pass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp10, gq15 pass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq15.pass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp10, gq20 pass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq20.pass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp10, gq20 NOpass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq20.no.filterpass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt




# echo $(cat vars_dropped_in_QC_zhaoming.txt| paste -s -d"|")
# wc -l vars_dropped_in_QC_zhaoming.txt
# 67 
# for var in $(cat vars_dropped_in_QC_zhaoming.txt); do
# 	echo "doing: $var"
# 	zgrep -w $var MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz| awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' 
# done
TYPE="zhaoming"
var="chr5:112839210:AG:A|chr17:43094316:TGATTCAGACTCCCCATCATGTGAGTCATCAGAACCTAACA:T|chr17:43124027:ACT:A|chr17:43070931:TTCTTCTGGGGTCAGGCCAG:T|chr13:32339062:CAG:C|chr13:32339287:AAAAG:A|chr13:32339482:TTATG:T|chr13:32337160:TAAAC:T|chr13:32337160:TA:T|chr16:68823523:CTG:C|chr14:95103391:G:GT|chr14:95117713:TA:T|chr2:47799609:CAAAG:C|chr17:31227281:G:GA|chr17:31259066:CA:C|chr17:31229953:CAA:C|chr17:31258501:AG:A|chr16:23635228:AC:A|chr9:37020778:C:CA|chr9:95485696:G:GTA|chr9:95458145:GC:G|chr13:48345089:AT:A|chr13:48349012:TA:T|chr13:48381280:AT:A|chr13:48364922:TA:T|chr13:48303990:G:GC|chr10:102627189:AG:A|chr10:102630077:CAA:C|chr17:7676152:C:CG|chr16:2060686:AC:A|chr11:32396281:G:GT|chr15:90765412:GATGTAAGTT:G|chr17:61715948:T:TA|chr17:8237439:GCTTT:G|chr11:10764306:AG:A|chr2:232210386:C:CT|chr19:45364832:CCTCA:C|chr2:127272934:CT:C|chr2:127286924:C:CACTGCT|chr2:127289348:CAA:C|chr13:102872217:C:CA|chr13:102862141:CAG:C|chr10:49482711:T:TC|chr10:49472935:CTGGGTCATAGA:C|chr10:49461381:CCT:C|chr16:89740841:GAGA:G|chr6:35457937:G:GC|chr11:22625192:CCGCTATCACCTTCAGGAAGTTGTTCTGAGGCAAG:C|chr9:35074172:GGACGGATCCA:G|chr14:45159189:C:CA|chr14:45189222:AC:A|chr14:45137187:GTGTT:G|chr4:54699759:ACAGT:A|chr5:132604019:CAA:C|chr5:132557417:T:TA|chr5:132609339:CAAAG:C|chr17:35101325:AGCCTCCC:A|chr8:144513108:CAT:C|chr8:144516695:GC:G|chr8:144514362:TC:T|chr16:3589547:CCT:C|chr8:31081140:AG:A|chr9:97689554:CA:C|chr3:14156338:GAGTAGACCGC:G|chr11:112093129:C:T|chr11:108325416:C:T|chr11:47234610:C:T"
echo -e "${TYPE} vars from VCF with no filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.vcf.gz| wc -l)"  >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp5, gq5 NOpass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq5.no.filter.pass.vcf.gz| wc -l)"  >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp5, gq5 pass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq5.pass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp5, gq10 pass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp5.gq10.pass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp10, gq15 pass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq15.pass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp10, gq20 pass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq20.pass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt
echo -e "${TYPE} vars from VCF after dp10, gq20 NOpass filter\t $(zgrep -Ew $var MERGED_biallelic_sorted_sjlife_1_2_${TYPE}_ID_edited.dp10.gq20.no.filterpass.vcf.gz | wc -l)" >> variant_counts_at_different_hard_filters.txt

# Most of these variants failed because they did not pass VQSR during hard filtering
# qin vars from VCF with no filter         108
# qin vars from VCF after dp5, gq5 NOpass filter   108
# qin vars from VCF after dp5, gq5 pass filter     1
# qin vars from VCF after dp5, gq10 pass filter    1
# qin vars from VCF after dp10, gq15 pass filter   1
# qin vars from VCF after dp10, gq20 pass filter   1
# qin vars from VCF after dp10, gq20 NOpass filter         108
# zhaoming vars from VCF with no filter    71
# zhaoming vars from VCF after dp5, gq5 NOpass filter      71
# zhaoming vars from VCF after dp5, gq5 pass filter        0
# zhaoming vars from VCF after dp5, gq10 pass filter       0
# zhaoming vars from VCF after dp10, gq15 pass filter      0
# zhaoming vars from VCF after dp10, gq20 pass filter      0
# zhaoming vars from VCF after dp10, gq20 NOpass filter    71
