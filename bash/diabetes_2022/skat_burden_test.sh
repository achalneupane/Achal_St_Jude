## Merge all plink files
# plink --bfile MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr1.PASS.decomposed_geno.0.1_hwe.1e-10 --make-bed --merge-list merge_all.list --keep-allele-order --out MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chrALL.PASS.decomposed_geno.0.1_hwe.1e-10


## Extract gnomAD with maf < 0.01 lines
head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt | sed 's/\t/\n/g' | nl
# head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01.txt

for i in {1..22}; do
CHROM="chr${i}"
echo "Doing ${CHROM}"
head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| awk '{print $178"\t"$179"\t"$181"\t"$182"\t"$38"\t"$39"\t"$44}' > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01_${CHROM}.txt	
awk -F'\t' '{if($38 != "." && $38 < 0.01) print $178"\t"$179"\t"$181"\t"$182"\t"$38"\t"$39"\t"$44}' "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHROM}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt" >> /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01_${CHROM}.txt
done

##############################################################################################
## Extract and merge all plink files with variants that are maf < 0.01 in gnomad_genome_ALL ##
##############################################################################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR

## EUR
for i in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${i}.PASS.decomposed_geno.0.1_hwe.1e-10 --extract extract_variants_gnomAD_ALL_NFE_lt_0.01_vars.list --keep EUR_samples.list --max-maf 0.01 --keep-allele-order --make-bed --out tmp_NFE_filtered_variants_EUR_chr${i}.dat
done

## AFR
for i in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${i}.PASS.decomposed_geno.0.1_hwe.1e-10 --extract extract_variants_gnomAD_ALL_AFR_lt_0.01_vars.list --keep AFR_samples.list --max-maf 0.01 --keep-allele-order --make-bed --out tmp_AFR_filtered_variants_AFR_chr${i}.dat
done

## Merge all EUR chromosomes
ls *EUR_chr*.dat.bim| sort -V| grep -v chr1.dat| sed 's/.bim//g' > merge_list_EUR.list
plink --bfile tmp_NFE_filtered_variants_EUR_chr1.dat --make-bed --merge-list merge_list_EUR.list --keep-allele-order --out tmp_NFE_filtered_variants_EUR_chr_ALL.dat
# Performing single-pass merge (3113 people, 21272 variants).
# Merged fileset written to tmp_NFE_filtered_variants_EUR_chr_ALL.dat-merge.bed +
# tmp_NFE_filtered_variants_EUR_chr_ALL.dat-merge.bim +
# tmp_NFE_filtered_variants_EUR_chr_ALL.dat-merge.fam .
# 21272 variants loaded from .bim file.
# 3113 people (0 males, 0 females, 3113 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to tmp_NFE_filtered_variants_EUR_chr_ALL.dat.nosex .
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 3113 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.998386.
# 21272 variants and 3113 people pass filters and QC.




## Merge all AFR chromosomes
ls *AFR_chr*.dat.bim| sort -V| grep -v chr1.dat| sed 's/.bim//g' > merge_list_AFR.list
plink --bfile tmp_AFR_filtered_variants_AFR_chr1.dat --make-bed --merge-list merge_list_AFR.list --keep-allele-order --out tmp_AFR_filtered_variants_AFR_chr_ALL.dat
# Performing single-pass merge (575 people, 26910 variants).
# Merged fileset written to tmp_EUR_filtered_variants_AFR_chr_ALL.dat-merge.bed +
# tmp_EUR_filtered_variants_AFR_chr_ALL.dat-merge.bim +
# tmp_EUR_filtered_variants_AFR_chr_ALL.dat-merge.fam .
# 26910 variants loaded from .bim file.
# 575 people (0 males, 0 females, 575 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to tmp_EUR_filtered_variants_AFR_chr_ALL.dat.nosex .
# Using 1 thread (no multithreaded calculations invoked.
# Before main variant filters, 575 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.998034.
# 26910 variants and 575 people pass filters and QC.



## Test whether --freq counts gives corrects counts for MAC
 # head  EUR_samples.list  > test_EUR_samples.list
 head extract_variants_gnomAD_ALL_NFE_lt_0.01_vars.list > test_EUR_vars.list

plink --bfile tmp_NFE_filtered_variants_EUR_chr_ALL.dat --extract test_EUR_vars.list --keep-allele-order --recodeA --out test_dat
plink --bfile tmp_NFE_filtered_variants_EUR_chr_ALL.dat --extract test_EUR_vars.list --keep-allele-order --make-bed --out test_dat_bed
plink2 --bfile test_dat_bed --freq counts --out test_MAC_EUR
# Based on this above test data, the MAC was correctly calculated



###########################################
## calculate MAC in EUR and AFR datasets ##
###########################################
module load plink2
plink2 --bfile tmp_NFE_filtered_variants_EUR_chr_ALL.dat --freq counts --out MAC_EUR
plink2 --bfile tmp_AFR_filtered_variants_AFR_chr_ALL.dat --freq counts --out MAC_AFR

## R code extract_plink_subset.r step 2...

# Now extract 
# Extract final variants for gene-based analysis 
plink --bfile tmp_NFE_filtered_variants_EUR_chr_ALL.dat --extract extract_variants_gnomAD_ALL_NFE_lt_0.01_MAC_gt_3_at_least_2_vars_per_gene.list --keep-allele-order --make-bed --out EUR_diabetes_chrALL.dat 

plink --bfile tmp_AFR_filtered_variants_AFR_chr_ALL.dat --extract extract_variants_gnomAD_ALL_AFR_lt_0.01_MAC_gt_3_at_least_2_vars_per_gene.list --keep-allele-order --make-bed --out AFR_diabetes_chrALL.dat 



## These are the two gene sets for AFR and EUR gene-based analysis. These were P/LP+LoF variants with maf  < 0.01 in gnomad all and Ethnicity specific gnomad (EUR = gnomad_gnome_NFE; AFR= gnomad_genom_AFR), MAC of 3 or more, and also with at least two variants in each gene
extract-variants-final-EUR.GENE-SNP
extract-variants-final-AFR.GENE-SNP
# CHR GENE SNPID


##############
## EUR test ##
##############
### Now that we have the generated BFILE and RAW file, we have to create diretories per each CHR and place in each directory a RAW file and generate gene-setrs within each CHR dir
mkdir gene-based-analysis
cd gene-based-analysis
mkdir EUR
cd EUR


ln -s ../../EUR_diabetes_chrALL.dat.* .
ln -s ../../extract-variants-final-EUR.GENE-SNP .
# 1- CREATE DIRECTORY FOR EACH CHR
for chr in {1..22}; do  mkdir chr${chr}; done


# 2 - MAKE RAW FILES FOR EACH CHR AND MOVE THEM INTO EACH DIRECTORY
CHR_TO_PROCESS="$(seq 1 22)"
THREADS=4
module load parallel
# BFILE="EOAD-2828-logistic-BETA"
BFILE="EUR_diabetes_chrALL.dat"
parallel -j ${THREADS} plink --bfile ${BFILE} --keep-allele-order --allow-no-sex --recode A --output-chr MT --chr {} --out chr{}/${BFILE}-chr{} ::: ${CHR_TO_PROCESS}

# for chr in {1..22}; do
# /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/plink --bfile ${BFILE} --keep-allele-order --allow-no-sex --recode A --output-chr MT --chr ${chr} --out chr${chr}/${BFILE}-chr${chr}
# done


# 3 - Extract annot info per chr and create gene-subsets
for file in $(ls *.GENE-SNP); do
INFILE=${file}
for chr in chr{1..22}; do
awk -v chr=${chr} '{if($1==chr) print $2,$3}' ${INFILE} > ${chr}/${chr}-${INFILE}
done
done



# 4 - Generate gen-sets using script get_gene.perl
set="EUR"    

BASE="extract-variants-final"
HOMEDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/EUR/"
for chr in {1..22}; do \
cd chr${chr}; mkdir geneset-${set}; perl ../../get_gene.perl chr${chr}-${BASE}-${set}.GENE-SNP; mv *.gene geneset-${set}; cd ${HOMEDIR}; done

# Now Run SKAT-loop-OPTIMAL.r

##############
## AFR test ##
##############
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis
mkdir AFR
cd AFR


ln -s ../../AFR_diabetes_chrALL.dat.* .
ln -s ../../extract-variants-final-AFR.GENE-SNP .
# 1- CREATE DIRECTORY FOR EACH CHR
for chr in {1..22}; do  mkdir chr${chr}; done


# 2 - MAKE RAW FILES FOR EACH CHR AND MOVE THEM INTO EACH DIRECTORY
CHR_TO_PROCESS="$(seq 1 22)"
THREADS=4
module load parallel
# BFILE="EOAD-2828-logistic-BETA"
BFILE="AFR_diabetes_chrALL.dat"
parallel -j ${THREADS} plink --bfile ${BFILE} --keep-allele-order --allow-no-sex --recode A --output-chr MT --chr {} --out chr{}/${BFILE}-chr{} ::: ${CHR_TO_PROCESS}


# for chr in {1..22}; do
# /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/plink --bfile ${BFILE} --keep-allele-order --allow-no-sex --recode A --output-chr MT --chr ${chr} --out chr${chr}/${BFILE}-chr${chr}
# done




# 3 - Extract annot info per chr and create gene-subsets
for file in $(ls *.GENE-SNP); do
INFILE=${file}
for chr in chr{1..22}; do
awk -v chr=${chr} '{if($1==chr) print $2,$3}' ${INFILE} > ${chr}/${chr}-${INFILE}
done
done


# 4 - Generate gen-sets using script get_gene.perl
set="AFR"    

BASE="extract-variants-final"
HOMEDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/AFR/"
for chr in {1..22}; do \
cd chr${chr}; mkdir geneset-${set}; perl ../../get_gene.perl chr${chr}-${BASE}-${set}.GENE-SNP; mv *.gene geneset-${set}; cd ${HOMEDIR}; done


# Now Run SKAT-loop-OPTIMAL.r




## gnomAD allele frequency for variants from Cindy (Email: 12/3/2022 2:58 pm)
DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/"
head -1 ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > gnomAD_lookup_list_Cindy_12_06_2022_out.txt
IFS=$'\n' # set IFS
for line in $(cat gnomAD_lookup_list_Cindy_12_06_2022.txt); do
CHR="$(echo ${line}| cut -d$'\t' -f1)"
BP="$(echo ${line}| cut -d$'\t' -f2)"
echo $CHR
grep -w $BP ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt >> gnomAD_lookup_list_Cindy_12_06_2022_out.txt
done




# Cindy's/ Yadav's email: 1/10/2023
# Can you please do the following for the diabetes GWAS? I have copied GWAS results for EUR (chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA) and AFR (chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA) in /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes.

## Run Mr-Mega
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/mr_mega/MR-MEGA -h

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega/v2_with_preQC_VCF


R --slave --vanilla --args input=chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA out=outputfilename < /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/mr_mega/fixP.r


PHENO: /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/SJLIFE_T2D_GWAS_EUR.pheno
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/SJLIFE_T2D_GWAS_AFR.pheno


ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf.gz.* .


## Calculate MAF
for CHR in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz --keep AFR_samples_v1.list --freq --out tmp_AFR_filtered_variants_AFR_freq_chr${CHR}
done



for CHR in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz --keep EUR_samples_v1.list --freq --out tmp_EUR_filtered_variants_EUR_freq_chr${CHR}
done


# Run Mr-Mega
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/mr_mega/MR-MEGA

## RUn on all variants in AFR and EUR
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/mr_mega/MR-MEGA -i MR-MEGA.in --gc --pc -1 -o MEGA.FINAL.RESULTS.txt
## Run on overlapping variants only
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/mr_mega/MR-MEGA -i mrmega.input --gc --pc -1 -o MEGA.FINAL.RESULTS_overlapping_vars.txt


# ## Request 2 
# No problem, thanks for your help! To clarify:

# (1) Get all variants with P<5e-6 in the AFR-only analysis
# (2) Match up variants in (1) with variants also present in the EUR analysis (i.e., replication): put these in an excel file
# (3) For all variants in (1), please send a rdata file with individual-level data for both AFR and EUR SJLIFE participants (I will find the email I got from Yadav and re-forward it to you so you can match up the variables included)

# Please let me know if this makes sense.

# Thanks,
# Cindy


# for CHR in {1..22}; do
# awk 'NR==FNR{a[$0];next} {for (i in a) if ($0 ~ i) {print}}' TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS.txt /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz.bim >>  TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS_BIM_extract.txt
# done

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega/v2_with_preQC_VCF

cat TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS.txt| parallel -j25 'chr=$(grep {} TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS.txt|sed "s/chr\([0-9]*\).*/\1/"); cat /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.preQC_biallelic_renamed_ID_edited.vcf.gz.bim | grep -w {}' >> bim_matched_TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS.txt
cut -d$'\t' -f2 bim_matched_TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS.txt > extract_snp.list

for chr in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.preQC_biallelic_renamed_ID_edited.vcf.gz --extract extract_snp.list --make-bed --out extract_snp_chr${chr}
done


plink --bfile extract_snp_chr1 --merge-list merge_list --keep-allele-order --out PLINK_TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS

## Check which variants match with R

## Recode 
plink --bfile PLINK_TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS --keep-allele-order --recode A --out PLINK_TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS_recodeA



## Try Metasoft
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/METASOFT
module load R
module load java/1.8.0_66
./plink2metasoft.py gwas2metasoft_file chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA
java -jar /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/METASOFT/software-notes-master/docs/files/METASOFT/Metasoft.jar --help

## Use R code: Chunk 4 to create input file
## METASOFT analysis
java -jar /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/METASOFT/software-notes-master/docs/files/METASOFT/Metasoft.jar -input METASOFT_input_8899222_overlapping_variants.meta -output metasoft_res
sed '1s/\#//g' metasoft_res > metasoft_res_edited


## Extract 185 SNP not previously sent Rdata, but in meta-analysis with P<5e-06
for chr in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.preQC_biallelic_renamed_ID_edited.vcf.gz --extract extract_SNPs.txt --make-bed --out extract_snp_chr${chr}
done


plink --bfile extract_snp_chr1 --merge-list merge_list --keep-allele-order --out PLINK_5e-06_meta-analysis_not_in_previous_RDATA

## Check which variants match with R

## Recode 
plink --bfile PLINK_5e-06_meta-analysis_not_in_previous_RDATA --keep-allele-order --recode A --out PLINK_5e-06_meta-analysis_not_in_previous_RDATA_recodeA

## Folloup_jan_24_2023
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/followup_jan_24_2023

for chr in 5 12; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.preQC_biallelic_renamed_ID_edited.vcf.gz --extract extract_SNP.list --make-bed --out extract_snp_chr${chr}
done

plink --bfile extract_snp_chr5 --bmerge extract_snp_chr12 --keep-allele-order --out followup_24_2023



plink --bfile followup_24_2023 --keep AFR_samples.list --freq --out followup_24_2023_AFR_freq
plink --bfile followup_24_2023 --keep EUR_samples.list --freq --out followup_24_2023_EUR_freq



## Get individual genotype data
plink --bfile followup_24_2023 --keep AFR_samples.list --keep-allele-order --make-bed --out folowup_24_2023_AFR
plink --bfile followup_24_2023 --keep EUR_samples.list --keep-allele-order --make-bed --out folowup_24_2023_EUR

plink --bfile folowup_24_2023_AFR --keep-allele-order --recode A --out folowup_24_2023_AFR_recodeA
plink --bfile folowup_24_2023_EUR --keep-allele-order --recode A --out folowup_24_2023_EUR_recodeA


## Locus Zoom
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/followup_jan_24_2023/locuszoom
ln -s ../../chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA .
ln -s ../../chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA .

# grep chr5 only and get : chr5_PCA_eur.assoc.logistic.clean.Psorted.withBETA.locus and chr5_PCA_afr.assoc.logistic.clean.Psorted.withBETA.locus

# with open("chr5_PCA_eur.assoc.logistic.clean.Psorted.withBETA.locus", "r") as f:
#     lines = f.readlines()
#     with open("chr5_PCA_eur.assoc.logistic.clean.Psorted.withBETA.locus.txt", "w") as out:
#         for line in lines:
#             elements = line.strip().split(" ")
#             new_line = "chr" + elements[0] + "\t" + elements[1] + "\t" + "chr" + elements[0] + ":" + elements[1] + "\t" + elements[13] + "\t" + elements[3] + "\t" + elements[11] + "\n"
#             out.write(new_line)

awk '{printf "chr%s\t%s\tchr%s:%s\t%s\t%s\t%s\n", $1, $3, $1, $3, $14, $4, $12}'  chr5_PCA_afr.assoc.logistic.clean.Psorted.withBETA.locus > chr5_PCA_afr.assoc.logistic.clean.Psorted.withBETA.locus.txt
awk '{printf "chr%s\t%s\tchr%s:%s\t%s\t%s\t%s\n", $1, $3, $1, $3, $14, $4, $12}'  chr5_PCA_eur.assoc.logistic.clean.Psorted.withBETA.locus > chr5_PCA_eur.assoc.logistic.clean.Psorted.withBETA.locus.txt

awk 'BEGIN {FS=OFS="\t"} NR==1 {print; next} {print $0 | "sort -k2,2n"}' chr5_PCA_afr.assoc.logistic.clean.Psorted.withBETA.locus.txt> locus_zoom_chr5_afr
awk 'BEGIN {FS=OFS="\t"} NR==1 {print; next} {print $0 | "sort -k2,2n"}' chr5_PCA_eur.assoc.logistic.clean.Psorted.withBETA.locus.txt> locus_zoom_chr5_eur
