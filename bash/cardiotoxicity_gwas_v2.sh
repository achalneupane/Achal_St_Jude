#!/usr/bin/bash


################
## Analysis 2 ##
################
# Yadav:  I think there are some SNPs which are probably <1% in controls but higher in cases and this the overall freq is higher than 1%. Can you remove variants from the summary statistics file that have <1% in controls only and then condition for chr16:25611595:C:T only?

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2
ln -s ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz .
PHENO=../pheno/sjlife_eur_dox_only_pcs.pheno

# extract CMP samples from $PHENO and variants from ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb

awk '{print $1"\t"$2}' $PHENO > samples.list

module load plink/1.90b
plink --chr 16 --from-bp 25361595 --make-bed --out sjlife_CMP --to-bp 25861595 --vcf MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz


awk '{print $2}' ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb > samplesnp.list
# keep samples and variants
ln -s ../finemap/sjlife_CMP.* .
plink --bfile sjlife_CMP --extract samplesnp.list --keep samples.list --keep-allele-order --make-bed --out samplesnp.dat


# plink --bfile sjlife_CMP --freq --keep-allele-order --out sjlife.freq.out
# awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' sjlife.freq.out.frq > sjlife.freq.out.frq_edited1

plink --bfile samplesnp.dat --freq --keep-allele-order --out samplesnp.dat.freq.out
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' samplesnp.dat.freq.out.frq > samplesnp.dat.freq.out.frq_edited1


# Yadav:  I think there are some SNPs which are probably <1% in controls but higher in cases and this the overall freq is higher than 1%. Can you remove variants from the summary statistics file that have <1% in controls only and then condition for chr16:25611595:C:T only?
plink --bfile  samplesnp.dat --keep controls.list --keep-allele-order --make-bed --out samplesnp.dat.controls
plink --bfile samplesnp.dat.controls --freq --keep-allele-order --out samplesnp.dat.controls.freq.out
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' samplesnp.dat.controls.freq.out.frq > samplesnp.dat.controls.freq.out.frq_edited1


plink --bfile sjlife_CMP --extract final.snp.list --keep samples.list --keep-allele-order --make-bed --out samplesnp.dat.final

## Create LD matrix
plink --r2 --bfile samplesnp.dat --matrix --out samplesnp_ld_matrix


# # chmod 755 *
# module load dos2unix
# # dos2unix samplesnp.z

##################################
## Run Finemap genomic region 1 ##
##################################
/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 --sss --log --in-files samplesnp --dataset 1

######################
## Run GEC analysis ##
######################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2
module load java
java -jar -Xmx1g /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/gec/gec.jar --var-gec --effect-number --plink-binary samplesnp.dat --genome --out samplesnp_chr16_vars.dat_GEC_test1
# java -jar -Xmx1g gec.jar --var-gec --effect-number --plink-binary samplesnp.dat --maf 0.01 --genome --out samplesnp_chr16_vars.dat_GEC_test1
# calculate Bonferroni corrected P-value based on Effective_Number in samplesnp_chr16_vars.dat_GEC_test1.sum
Bonferroni_P_chr16 = 0.0002899896 = 0.05/172.42 

###########################
## run GCTA analysis SNP ##
###########################
module load gcta
gcta64  --bfile samplesnp.dat --chr 16 --cojo-file samplesnp_vars.ma --cojo-slct --cojo-p 0.0002899896 --cojo-wind 100000 --out chr16_cojo
# gcta64  --bfile samplesnp.dat --chr 16 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.00017925 --cojo-wind 100000 --out TITN_cojo
# Finally, 4 associated SNPs are selected.

# Running cojo with selected SNP list 
gcta64  --bfile samplesnp.dat --chr 16 --cojo-file samplesnp_vars.ma --cojo-p 0.0002899896 --cojo-cond chr16_snplist --out chr16_snplist_coho.out
# gcta64  --bfile samplesnp.dat --chr 16 --maf 0.01 --cojo-file samplesnp_vars.ma --cojo-p 0.0002899896 --cojo-cond chr16_snplist --out chr16_25611595_snplist.out

# # <chr16_snplist>
# chr16:25554098:G:A
# chr16:25573423:G:GTA
# chr16:25609966:C:T
# chr16:25611595:C:T

# Yadav : "Can you change cojo-p to 5e-08?"
gcta64  --bfile samplesnp.dat --chr 16 --cojo-file samplesnp_vars.ma --cojo-slct --cojo-p 5e-08 --cojo-wind 100000 --out chr16_cojo_5e_08
# Finally, 2 associated SNPs are selected...

# chr16:25573423:G:GTA
# chr16:25609966:C:T


# Yadav: "Can you just condition on chr16:25611595:C:T and see if there are additional SNPs with P<5e-08?"
gcta64  --bfile samplesnp.dat --chr 16 --cojo-file samplesnp_vars.ma --cojo-p 5e-08 --cojo-cond chr_snplist_chr16_25611595 --out chr16_snplist_chr16_25611595_coho.out
# result in spreadsheet tab: cojo_conditioned_chr16_25611595

###################################################
## Run Finemap with assumption for 4 causal SNPs ##
###################################################
mkdir finemap_with_4_associations
cd finemap_with_4_associations
ln -s ../samplesnp.z .
ln -s ../samplesnp .
ln -s ../samplesnp_ld_matrix.ld .
/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 --config --in-files samplesnp --dataset 1 --rsids chr16:25554098:G:A,chr16:25573423:G:GTA,chr16:25609966:C:T,chr16:25611595:C:T
# |--------------------------------------|
# | Welcome to FINEMAP v1.4.1            |
# |                                      |
# | (c) 2015-2022 University of Helsinki |
# |                                      |
# | Help :                               |
# | - ./finemap --help                   |
# | - www.finemap.me                     |
# | - www.christianbenner.com            |
# |                                      |
# | Contact :                            |
# | - finemap@christianbenner.com        |
# | - matti.pirinen@helsinki.fi          |
# |--------------------------------------|

# --------
# SETTINGS
# --------
# - dataset     : 1
# - corr-config : 0.95
# - prior-std   : 0.05

# --------------------
# CAUSAL CONFIGURATION
# --------------------
# - GWAS summary stats               : samplesnp.z
# - SNP correlations                 : samplesnp_ld_matrix.ld
# - Reading input                    : done!

# - Updated prior SD of effect sizes : 0.05 0.0963 0.185 0.357

# - SNPs in causal configuration     : chr16:25554098:G:A,chr16:25573423:G:GTA,chr16:25609966:C:T,chr16:25611595:C:T
# ------------------------------------
# - Log-Likelihood of configuration  : -1085.397951
# - ML causal effect estimates       : 0.398051,0.480682,1.69795,0.490914
# - SE of ML causal effect estimates : 0.0709777,0.0652414,0.0831256,0.0630158
# - P-values for causal effects      : 2.04545e-08,1.73546e-13,9.75993e-93,6.68381e-15
# ------------------------------------
# - Log10-BF of configuration        : 133.829703
# - Posterior mean of causal effects : 0.37346,0.451057,1.59116,0.494195
# - Posterior SD of causal effects   : 0.0687502,0.0631939,0.0793122,0.060125
# ------------------------------------
# - Regional SNP heritability        : 0.435 (SD: 0.0274 ; 95% CI: [0.384,0.487])

# - Run time                         : 0 hours, 0 minutes, 0 seconds


