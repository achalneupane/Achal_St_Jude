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

