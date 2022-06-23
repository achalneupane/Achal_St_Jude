#!/bin/bash
### ##################
## Date: 06/06/2022 ##
## Achal Neupane    ##
######################
#################################
## Phenotype data manipulation ##
#################################
# Wrote an R code : phenotype_cleaning_st_jude_life.r

##########################################################
## Script to merge VCF files from disjoint sample lists ##
##########################################################
# Count samples in VCF
# zcat $VCF |head -5000|  grep "#CHROM" | tr "\t" "\n " | tail -n +10 | uniq | wc -l

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC

ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife_1/preQC/* .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife_2/preQC/* .
(module avail) |& grep -i bcftools| sort -V
chmod 755 entrypoint_merge_vcf.sh
dos2unix entrypoint_merge_vcf.sh

for CHR in {1..22}; do \
	unset VCF1; unset VCF2; unset MERGED; \
	echo "Doing chr${CHR}"; \
	export VCF1="SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${CHR}.vcf.gz"; \
	export VCF2="SJLIFE2.GATKv3.4.VQSR.sjlid_chr${CHR}.vcf.gz"; \
	echo -e "**\nMerging ${VCF1} \nand \n${VCF2}\n**"; \
	export MEM=6; \
	export THREADS=4; \
	export MERGED="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC"; \
	export OUT_DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/"; \
	bsub \
	-P "chr${CHR}_merge" \
	-J "chr${CHR}_merge" \
	-o "${OUT_DIR}/logs/${MERGED}_s00VCFmerge.%J" \
	-n ${THREADS} \
	-R "rusage[mem=8192]" \
	"./entrypoint_merge_vcf.sh"; \
done; 



## Convert to biallelic; reset IDs; reset sample names
mkdir -p /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned
ln -s ../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC.vcf.gz .
cat ../rename_VCF_ids_to_sjlids.txt ../SJLIFEWGS_3006_vcfid_to_sjlid.linkfile > rename_link_file.txt .
mkdir logs

##########################
##########################
##########################

# -f (reference sequence) in bcftools will turn on left-alignment and normalization, however, see also the --do-not-normalize option below. (https://samtools.github.io/bcftools/bcftools.html)
cat <<\EoF > edit_vcf_header.sh
#!/usr/bin/bash
module load bcftools/1.10.2

cd "${OUT_DIR}"
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz -Oz -o MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic.vcf.gz
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic.vcf.gz -Oz -o MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_ID_edited.vcf.gz
rm MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic.vcf.gz
bcftools reheader -s rename_link_file.txt MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_ID_edited.vcf.gz > MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz
rm MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_ID_edited.vcf.gz
bcftools index -f -t --threads 4 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz
EoF

for i in {1..22}; do \
	unset CHR; \
	export CHR="$i"; \
	echo "Doing chr${CHR}"; \
	export MEM=6; \
	export THREADS=4; \
	export OUT_DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/"; \
	bsub \
	-P "chr${CHR}_edit" \
	-J "chr${CHR}_edit" \
	-o "${OUT_DIR}/logs/chr${CHR}_edited.%J" \
	-n ${THREADS} \
	-R "rusage[mem=8192]" \
	"./edit_vcf_header.sh"; \
done; 

# #############################
# ## VCF to plink conversion ##
# #############################
# for CHR in {1..22}; do \
# 	unset VCF;unset PLINK_FILE; \
# 	echo "Doing chr${CHR}"; \
# 	export VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz"; \
# 	export PLINK_FILE="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${CHR}.PASS.decomposed"; \
# 	export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/"
# 	echo -e "**\nConverting ${VCF} \n to Plink \n${PLINK_FILE}\n**"; \
# 	export MEM=6; \
# 	export THREADS=4; \
# 	bsub \
# 	-P "chr${CHR}_plink" \
# 	-J "chr${CHR}_plink" \
# 	-o "${WORKDIR}/logs/${PLINK_FILE}_s01VCFplink.%J" \
# 	-n ${THREADS} \
# 	-R "rusage[mem=8192]" \
# 	"./entrypoint_VCF_to_Plink_with_basic_QC.sh"; \
# done; 





@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
####################
## VCF annotation ##
####################
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
ln -s ../../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz* .
## See VCF_merge_annotate_QC.sh for details about generating reference files and databases

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
	-P "${CHR}_annotate" \
	-J "${CHR}_annotate" \
	-o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
	-n ${THREADS} \
	-R "rusage[mem=8192]" \
	"./entrypoint_VCFannotation.sh"; \
done; 



# #####################################
# ## Annovar and intervar Annotation ##
# #####################################
# # First check avblist hg38_avdblist.txt for the latest releases
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar
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


# ## Annovar annotation ##
# # https://annovar.openbioinformatics.org/en/latest/user-guide/startup/
# # NOTE: the -operation argument tells ANNOVAR which operations to use for each of the protocols: g means gene-based, gx means gene-based with cross-reference annotation (from -xref argument), r means region-based and f means filter-based.

# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/ANNOVAR_ANNOTATION
# ln -s ../../*.vcf.gz .

# ## Intervar ##
# ## Add Intervar annotation
# # --input_type AVinput or VCF or VCF_m

# for i in {1..22}; do \
# export CHR="chr${i}"; \
# echo "Annotating $CHR"; \
# unset INPUT_VCF; \
# export THREADS=4; \
# export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.PASS.decomposed.vcf.gz"; \
# export JAVA="java"; \
# export JAVAOPTS="-Xms4g -Xmx30g"; \
# export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/ANNOVAR_ANNOTATION/"; \
# export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
# export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
# export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
# export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
# export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
# bsub \
# -P "${CHR}_ANNOVAR" \
# -J "${CHR}_ANNOVAR" \
# -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_ANNOVAR.%J" \
# -n ${THREADS} \
# -R "rusage[mem=8192]" \
# "./entrypoint_ANNOVAR_annotation.sh"; \
# done;


##########################
## Edit PRS score files ##
##########################

# Sites to look for PRS:
# https://prsweb.sph.umich.edu:8443/displayData/displayTable?select_desc=174.1&select_phenomes=MGI&select_odds=2
# https://www.pgscatalog.org/score/PGS000015/
##################################
## Breast cancer (Michigan web) ##
##################################
awk ' {FS=OFS="\t"; $(NF+1) = "Chr"$1":"$2 }1' PRSWEB_PHECODE174.1_Onco-iCOGS-Overall-BRCA_PRS-CS_UKB_20200608_WEIGHTS.txt | grep -v ^## > OVERALL_PRSWEB_PHECODE174.1_Onco-iCOGS-Overall-BRCA_PRS-CS_UKB_20200608_WEIGHTS_edited_1.txt
awk ' {FS=OFS="\t"; $(NF+1) = "Chr"$1":"$2 }1' PRSWEB_PHECODE174.1_Onco-iCOGS-ER-positive-BRCA_PRS-CS_MGI_20200608_WEIGHTS.txt | grep -v ^## > ER_POS_PRSWEB_PHECODE174.1_Onco-iCOGS-ER-positive-BRCA_PRS-CS_MGI_20200608_WEIGHTS_edited_1.txt
awk ' {FS=OFS="\t"; $(NF+1) = "Chr"$1":"$2 }1' PRSWEB_PHECODE174.1_GWAS-Catalog-r2019-05-03-X174.1_PT_UKB_20200608_WEIGHTS.txt| grep -v ^## > ER_NEG_PRSWEB_PHECODE174.1_GWAS-Catalog-r2019-05-03-X174.1_PT_UKB_20200608_WEIGHTS_edited_1.txt


awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' OVERALL_PRSWEB_PHECODE174.1_Onco-iCOGS-Overall-BRCA_PRS-CS_UKB_20200608_WEIGHTS_edited_1.txt | sed -n '1d;p' | head> test.bed

# Print last col: awk '{print $(NF)}'
awk '{print "chr"$1"\t"$2"\t"$2+1"\t"$6}' OVERALL_PRSWEB_PHECODE174.1_Onco-iCOGS-Overall-BRCA_PRS-CS_UKB_20200608_WEIGHTS_edited_1.txt | sed -n '1d;p' > Michigan_OVERALL.bed
awk '{print "chr"$1"\t"$2"\t"$2+1"\t"$6}' ER_POS_PRSWEB_PHECODE174.1_Onco-iCOGS-ER-positive-BRCA_PRS-CS_MGI_20200608_WEIGHTS_edited_1.txt | sed -n '1d;p' > Michigan_ER_POS.bed
awk '{print "chr"$1"\t"$2"\t"$2+1"\t"$9}' ER_NEG_PRSWEB_PHECODE174.1_GWAS-Catalog-r2019-05-03-X174.1_PT_UKB_20200608_WEIGHTS_edited_1.txt | sed -n '1d;p' > Michigan_ER_NEG.bed

/home/aneupane/liftover/liftOver Michigan_OVERALL.bed /home/aneupane/liftover/hg19ToHg38.over.chain Michigan_OVERALL_GrCh38.bed Michigan_OVERALL_unmapped.bed
# [aneupane@noderome164 all_downloads]$ /home/aneupane/liftover/liftOver Michigan_OVERALL.bed /home/aneupane/liftover/hg19ToHg38.over.chain Michigan_OVERALL_GrCh38.bed Michigan_OVERALL_unmapped.bed
# Reading liftover chains
# Mapping coordinates


## Excel Match =INDEX(Mavaddat_2019_GrCH38_converted!A:A,MATCH(Mavaddat_2019_Table_S7_GRCh37!A2,Mavaddat_2019_GrCH38_converted!D:D,FALSE))
         # Something you want from different spreadsheet,MATCH(Match what?, Match Where?, exact match?
/home/aneupane/liftover/liftOver Michigan_ER_POS.bed /home/aneupane/liftover/hg19ToHg38.over.chain Michigan_ER_POS_GrCh38.bed Michigan_ER_POS_unmapped.bed
/home/aneupane/liftover/liftOver Michigan_ER_NEG.bed /home/aneupane/liftover/hg19ToHg38.over.chain Michigan_ER_NEG_GrCh38.bed Michigan_ER_NEG_unmapped.bed


###################
## Mavaddat 2019 ##
###################
/home/aneupane/liftover/liftOver Mavaddat_2019_table_S7_GrCh37.bed /home/aneupane/liftover/hg19ToHg38.over.chain Mavaddat_2019_table_S7_GrCh38.bed Mavaddat_2019_table_S7_unmapped.bed

################
## Khera 2018 ##
################
/home/aneupane/liftover/liftOver Khera_2018_Hg19.bed /home/aneupane/liftover/hg19ToHg38.over.chain Khera_2018_Hg19_GrCh38.bed Khera_2018_Hg19_unmapped.bed

####################
## Vijayakrishnan ##
####################
/home/aneupane/liftover/liftOver ALL_Vijayakrishnan_GRCh37.bed /home/aneupane/liftover/hg19ToHg38.over.chain ALL_Vijayakrishnan_Hg19_GrCh38.bed ALL_Vijayakrishnan_Hg19_unmapped.bed

##############################
## Pleiotropy all Positions ##
##############################
/home/aneupane/liftover/liftOver GRCh37_all_pos.bed /home/aneupane/liftover/hg19ToHg38.over.chain GRCh37_all_pos_Hg19_GrCh38.bed GRCh37_all_pos_Hg19_unmapped.bed

########################
## Wang et al African ##
########################
/home/aneupane/liftover/liftOver Wang_African_GRCh37.bed /home/aneupane/liftover/hg19ToHg38.over.chain Wang_African_GrCh38.bed Wang_African_Hg19_unmapped.bed



awk 'FNR==NR{a[$4] = (a[$4]==""?"":a[$4] " ") $2 OFS $3 OFS $4; next}
    {print $4, ($4 in a ? a[$4] : 0)}' Michigan_ER_NEG_GrCh38.bed ER_NEG_PRSWEB_PHECODE174.1_GWAS-Catalog-r2019-05-03-X174.1_PT_UKB_20200608_WEIGHTS_edited_1.txt


## Extract plink subset for PRS
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/plink_data
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz* .

cat ../SNVs_PRS.bed ../Indels_PRS.bed > PRS_vars.bed


module load bcftools/1.9
module load plink/1.90b

for it in {1..22}; do
echo "Doing chr${it}"; \
export CHR="$it"; \
export THREADS=4; \
	bsub \
	-P "chr${CHR}_extract" \
	-J "chr${CHR}_extract" \
	-e "${PWD}/logs/chr${CHR}_extractV3_err.%J" \
	-o "${PWD}/logs/chr${CHR}_extractV3.%J" \
	-n ${THREADS} \
	-R "rusage[mem=50000]" \
	"./extract_variants_from_VCF_for_PRS.sh"; \
done


# bjobs| grep RUN| grep extract| cut -d$' ' -f1| xargs bkill


## Merge plink files
for CHR in {1..22}; do
echo "PRS_chr${CHR}_v2" >> merge_list.txt
done
plink --make-bed --merge-list merge_list.txt  --out sjlife_all_PRS



## run this one more time for the variants that I missed in the first run
ln -s ../extract_batch2_SNVs_PRS.bed .
for it in 1 2 3 4 5 6 7 8 9 10 11 12 14 16 18 19 22; do
echo "Doing chr${it}"; \
export CHR="$it"; \
export THREADS=4; \
	bsub \
	-P "chr${CHR}_extract" \
	-J "chr${CHR}_extract" \
	-e "${PWD}/logs/chr${CHR}_extractV3_err.%J" \
	-o "${PWD}/logs/chr${CHR}_extractV3.%J" \
	-n ${THREADS} \
	-R "rusage[mem=20000]" \
	"./extract_variants_from_VCF_for_PRS_v3.sh"; \
done


## also run for missing meningioma variants
for it in 10 11; do
echo "Doing chr${it}"; \
export CHR="$it"; \
export THREADS=4; \
	bsub \
	-P "chr${CHR}_extract" \
	-J "chr${CHR}_extract" \
	-e "${PWD}/logs/chr${CHR}_extractV3_err.%J" \
	-o "${PWD}/logs/chr${CHR}_extractV3.%J" \
	-n ${THREADS} \
	-R "rusage[mem=20000]" \
	"./meningioma_extract_variants_from_VCF_for_PRS.sh"; \
done

## Merge plink files
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 14 16 18 19 22; do
echo "PRS_chr${CHR}_v3" >> merge_list2.txt
done
plink --make-bed --merge-list merge_list2.txt  --out sjlife_all_PRS2

# merge sjlife_all_PRS and sjlife_all_PRS2

plink --bfile sjlife_all_PRS --bmerge sjlife_all_PRS2 --out sjlife_all_PRS_all


for CHR in 10 11; do
echo "PRS_chr${CHR}_v4" >> merge_list3.txt
done
plink --make-bed --merge-list merge_list3.txt  --out sjlife_all_PRS3

# merge sjlife_all_PRS and sjlife_all_PRS2

plink --bfile sjlife_all_PRS_all --bmerge sjlife_all_PRS3 --out sjlife_all_PRS_all_final

# # Variants with PRS scores from all cancer types
# [aneupane@noderome186 plink_data]$ wc -l PRS_all_cancers_vars.txt
# 1546921 PRS_all_cancers_vars.txt


# ## Extract from plink file
# # Step 1.
# plink --bfile sjlife_all_PRS_all --extract PRS_all_cancers_vars.txt --make-bed --keep-allele-order --out sjlife_found_set1
# # 1166787 variants loaded from .bim file.
# # 4507 people (0 males, 0 females, 4507 ambiguous) loaded from .fam.
# # Ambiguous sex IDs written to sjlife_found_set1.nosex .
# # --extract: 1119983 variants remaining.
# # Using 1 thread (no multithreaded calculations invoked.
# # Before main variant filters, 4507 founders and 0 nonfounders present.
# # Calculating allele frequencies... done.
# # Total genotyping rate is 0.999865.
# # 1119983 variants and 4507 people pass filters and QC.
# # Note: No phenotypes present.
# # --make-bed to sjlife_found_set1.bed + sjlife_found_set1.bim +




# awk '{ print $2 }' sjlife_found_set1.bim > set1_vars

# # Step 2.
# # Excluding variants that did match in step 1
# plink --bfile sjlife_all_PRS_all --exclude set1_vars --make-bed --keep-allele-order --out sjlife_not_found
# # 1166787 variants loaded from .bim file.
# # 4507 people (0 males, 0 females, 4507 ambiguous) loaded from .fam.
# # Ambiguous sex IDs written to sjlife_not_found.nosex .
# # --exclude: 46804 variants remaining.
# # Using 1 thread (no multithreaded calculations invoked.
# # Before main variant filters, 4507 founders and 0 nonfounders present.
# # Calculating allele frequencies... done.
# # Total genotyping rate is 0.997626.
# # 46804 variants and 4507 people pass filters and QC.
# # Note: No phenotypes present.
# # --make-bed to sjlife_not_found.bed + sjlife_not_found.bim +
# # sjlife_not_found.fam ... done.


# # Step 3.
# awk '{ print $2 }' sjlife_not_found.bim > flip_list
# # Now flipping variants that did not match in step 1 (flip list obtained from step 2)
# plink --bfile sjlife_all_PRS_all --flip flip_list --make-bed --keep-allele-order --out sjlife_all_PRS_flipped
# # 1166787 variants loaded from .bim file.
# # 4507 people (0 males, 0 females, 4507 ambiguous) loaded from .fam.
# # Ambiguous sex IDs written to sjlife_all_PRS_flipped.nosex .
# # --flip: 45599 SNPs flipped.
# # Warning: 24882 variants had at least one non-A/C/G/T allele name.
# # Using 1 thread (no multithreaded calculations invoked.
# # Before main variant filters, 4507 founders and 0 nonfounders present.
# # Calculating allele frequencies... done.
# # Total genotyping rate is 0.999776.
# # 1166787 variants and 4507 people pass filters and QC.
# # Note: No phenotypes present.
# # --make-bed to sjlife_all_PRS_flipped.bed + sjlife_all_PRS_flipped.bim +





# ## Now extract again from the flipped. Here flipping alleles did not add any new matches
# plink --bfile sjlife_all_PRS_flipped --extract PRS_all_cancers_vars.txt --make-bed --keep-allele-order --out sjlife_found_v2
# # 1166787 variants loaded from .bim file.
# # 4507 people (0 males, 0 females, 4507 ambiguous) loaded from .fam.
# # Ambiguous sex IDs written to sjlife_found_v2.nosex .
# # --extract: 1119983 variants remaining.
# # Using 1 thread (no multithreaded calculations invoked.
# # Before main variant filters, 4507 founders and 0 nonfounders present.
# # Calculating allele frequencies... done.
# # Total genotyping rate is 0.999865.
# # 1119983 variants and 4507 people pass filters and QC.
# # Note: No phenotypes present.
# # --make-bed to sjlife_found_v2.bed + sjlife_found_v2.bim + sjlife_found_v2.fam



awk 'NR==FNR{a[$1":"$2]=$3" "$4;next}($1":"$4 in a){print $0, a[$1":"$4]}' ../ALL_Cancers_PRS_data.txt sjlife_all_PRS_flipped.bim | awk '($5==$7 || $5==$8) && ($6==$7 || $6==$8)' | awk '!a[$1":"$4]++' | wc -l
1119983

Of the 2241696 variants in the PRS file (which seem to have duplicates as well), there are 1121122 unique variants.
wc -l ../ALL_Cancers_PRS_data.txt
2241697 ../ALL_Cancers_PRS_data.txt
awk '!a[$1":"$2]++' ../ALL_Cancers_PRS_data.txt | wc -l
1121123







## Extract Clinvar, MetaSVM and LoF variants from VCF
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/annotated_variants/annotated_vars_from_PreQC_VCF_Clinvar_MetaSVM_LoF
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz* .

for it in {1..22}; do
echo "Doing chr${it}"; \
export CHR="$it"; \
export THREADS=4; \
	bsub \
	-P "chr${CHR}_extract" \
	-J "chr${CHR}_extract" \
	-e "${PWD}/logs/chr${CHR}_extract_err.%J" \
	-o "${PWD}/logs/chr${CHR}_extract.%J" \
	-n ${THREADS} \
	-R "rusage[mem=60000]" \
	"./extract_variants_from_VCF_for_Clinvar_MetaSVM_LoF_PreQC.sh"; \
done


## Merge all plink files
for CHR in {1..22}; do
echo "Annotated_Pathogenic_PreQC_chr${CHR}" >> merge_list.txt
done
plink --make-bed --merge-list merge_list.txt  --out sjlife_all_pathogenic_variants_PreQC
## Recode
plink --bfile sjlife_all_pathogenic_variants_PreQC --recodeA --out sjlife_all_pathogenic_variants_PreQC_recodeA
plink --bfile sjlife_all_pathogenic_variants_PreQC --recode12 --out sjlife_all_pathogenic_variants_PreQC_recode12
plink --bfile sjlife_all_pathogenic_variants_PreQC --recodeAD --out sjlife_all_pathogenic_variants_PreQC_recodeAD

# https://zzz.bwh.harvard.edu/plink/dataman.shtml
# plink --file data --recodeAD
# which, assuming C is the minor allele, will recode genotypes as follows:
#  SNP       SNP_A ,  SNP_HET
#      ---       -----    -----
#      A A   ->    0   ,   0
#      A C   ->    1   ,   1
#      C C   ->    2   ,   0
#      0 0   ->   NA   ,  NA

# n otherwords, the default for the additive recoding is to count the number of minor alleles per person. The --recodeAD option produces both an additive and dominance coding: use --recodeA instead to skip the SNP_HET coding.
# The --recodeAD option saves the data to a single file

## Therefore, --recodeA is the right option to get the carriers
####################################################################
## check variants from  SNPeff in Zhaoming and Qin et al's papers ##
####################################################################
# For this task, I am mostly using R code variants_counts.R
# First create a list of variants to be searched

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/variants_in_VCF
# ln -s ../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf .

# module load parallel
# getVars()
# {
# VCF="$1"
# CHR="$2"
# grep -v "^##" ${VCF} | cut -f3 | cut -d";" -f1 > CHR${CHR}_Vars.vcf
# }

# export -f getVars
# parallel -j22 getVars MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr{}.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf {} ::: {1..22}


# Then run R variants_counts_preQC_VCF_06_13_2022.R

##################
## LoF variants ##
##################
# Qin et al: frameshift      missense      nonsense    proteinDel        splice splice_region 
# 196            70           175             1            74             4 

## Annotation of interest
# frameshift_variant, start_lost, stop_gained, splice, ^splice_region_variant$
head -1 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt > LoF_variants.txt
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
