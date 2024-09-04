#!/bin/bash

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3
# To calculate maf, use sjlife.fam on /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3. Use pheno/sjlife_ttn_bag3.pheno
awk '{print $1"\t"$2}' pheno/sjlife_ttn_bag3.pheno > samples_for_maf.txt
awk '{print $1"\t"$2}' pheno/ccss_org_eur_cardiotoxic_exposed.pheno > samples_for_maf_ccss_org.txt
awk '{print $1"\t"$2}' pheno/ccss_exp_eur_cardiotoxic_exposed.pheno > samples_for_maf_ccss_exp.txt

# USE R code TITN_BAG3.r
# Chunk 1.

# LD calculation with plink; https://zzz.bwh.harvard.edu/plink/ld.shtml; see, "Alternatively, it is possible to add the --matrix option, which creates a matrix of LD values rather than a list: in this case, all SNP pairs are calculated and reported, even for SNPs on different chromosomes"
# Also: https://www.biostars.org/p/268141/
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64
ln -s ../../samples_for_maf.txt


module load plink/1.90b
plink --bfile sjlife --extract samplesnp.list --keep samples_for_maf.txt --keep-allele-order --make-bed --out samplesnp.dat

plink --bfile samplesnp.dat --keep samples_for_maf.txt --freq --keep-allele-order --out sjlife.freq.out

# First removing trailing spaces in the file, then removing chr
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' sjlife.freq.out.frq | awk -F " " '{gsub(/chr/, "", $2); print}' OFS="\t" > sjlife.freq.out.frq_edited1

## samplesnp.list
# chr2:178802888:A:G
# chr2:178803079:C:T
# chr2:178804123:G:A
# ...
# NOTE: The order changes according to the chr and position order in the SNP list. So, the summary statistics needs to be updated to maintain the order as in plink subset. 

# test1.list
# chr2:178313079:G:A
# chr2:178313658:T:C
# chr2:178313672:C:A
# chr2:178313675:A:G
# chr2:178313779:A:G
# chr2:178313920:A:G
# chr2:178314290:C:T
# chr2:178314441:T:C


# test2.list
# chr2:178313079:G:A
# chr2:178313658:T:C
# chr2:178314441:T:C
# chr2:178313920:A:G
# chr2:178313672:C:A
# chr2:178313675:A:G
# chr2:178313779:A:G
# chr2:178314290:C:T




# # ## Master file: <samplesnp>
# # z;ld;snp;config;cred;log;k;n_samples
# # samplesnp.z;samplesnp.ld;samplesnp.snp;samplesnp.config;samplesnp.cred;samplesnp_finemap.log;samplesnp.k;1645
# z;ld;snp;config;cred;log;n_samples
# samplesnp.z;samplesnp.ld;samplesnp.snp;samplesnp.config;samplesnp.cred;samplesnp_finemap.log;1645



# # I was getting errors with the file format, so I had to run dos2unix
# # chmod 777 *
# # dos2unix samplesnp.z

# ./finemap_v1.4.1_x86_64 --sss --log --in-files samplesnp --dataset 1
# # ./finemap_v1.4.1_x86_64 --sss --log --n-causal-snps 2 --in-files samplesnp
# ./finemap_v1.4.1_x86_64 --sss --log --in-files samplesnp



# # Additionally, .z file shoul have maf in specific form. I was getting this error: "Error : Expected a floating point value greater than 0.0 and smaller or equal than 0.5 in line 3 column 6 of file 'samplesnp.z'!"

## Extract TITN and BAG3 variants with MAF > 1%; separately.
plink --bfile samplesnp.dat --extract samplesnp_TITN_gt_MAF_1_perc_vars.list --keep-allele-order --make-bed --out samplesnp_TITN_gt_MAF_1_perc_vars.dat
plink --bfile samplesnp.dat --extract samplesnp_BAG3_gt_MAF_1_perc_vats.list --keep-allele-order --make-bed --out samplesnp_BAG3_gt_MAF_1_perc_vars.dat

# RUN GEC analysis to calculates the number of independent variants in datasets prepared above <http://pmglab.top/gec/#/download>; wget http://pmglab.top/gec/data/archive/v0.2/gecV0.2.zip
module load java
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/gec
ln -s ../samplesnp_TITN_gt_MAF_1_perc_vars.dat* .
ln -s ../samplesnp_BAG3_gt_MAF_1_perc_vars.dat* .

java -jar -Xmx1g gec.jar --var-gec --effect-number --plink-binary samplesnp_TITN_gt_MAF_1_perc_vars.dat --maf 0.01  --genome --out samplesnp_TITN_gt_MAF_1_perc_vars.dat_GEC_test1
java -jar -Xmx1g gec.jar --var-gec --effect-number --plink-binary samplesnp_BAG3_gt_MAF_1_perc_vars.dat --maf 0.01 --genome --out samplesnp_BAG3_gt_MAF_1_perc_vars.dat_GEC_test1
# calculate Bonferroni corrected P-value
BF_P_TITN = 0.00017925 = 0.05/278.94 
BF_P_BAG3 = 0.000103976 = 0.05/480.88

## Run conditional joint analysis; https://gcta.freeforums.net/thread/178/conditional-joint-analysis-using-summary
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test
ln -s ../samplesnp_TITN_gt_MAF_1_perc_vars.dat* .
ln -s ../samplesnp_BAG3_gt_MAF_1_perc_vars.dat* .

module load gcta
gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat  --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.00017925 --cojo-wind 100000 --out TITN_cojo
gcta64  --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat  --chr 10 --maf 0.01 --cojo-file samplesnp_BAG3_gt_MAF_1_perc_vats.ma --cojo-slct --cojo-p 0.000103976 --cojo-wind 100000 --out BAG3_cojo

## Test if GCTA options are valid with p 0.01 or even 0.05
# gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat  --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.05 --cojo-wind 100000 --out TITN_cojo
# gcta64  --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat  --chr 10 --maf 0.01 --cojo-file samplesnp_BAG3_gt_MAF_1_perc_vats.ma --cojo-slct --cojo-p 0.05 --cojo-wind 100000 --out BAG3_cojo

## Since gcta analysis did not find any independent SNPs, which implies we can assume 1 causal variant within each of the TTN and BAG3 locus

## Now running FINEMAP for each TITN and BAG3 assuming 1 causal variant
##########
## TITN ##
##########
# Create LD file
plink --r2 --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat --matrix --out samplesnp_TITN_gt_MAF_1_perc
plink --r2 --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat --matrix --out samplesnp_BAG3_gt_MAF_1_perc

## <samplesnp_TITN_gt_MAF_1_perc>
# z;ld;snp;config;cred;log;n_samples
# samplesnp_TITN_gt_MAF_1_perc.z;samplesnp_TITN_gt_MAF_1_perc.ld;samplesnp_TITN_gt_MAF_1_perc.snp;samplesnp_TITN_gt_MAF_1_perc.config;samplesnp_TITN_gt_MAF_1_perc.cred;samplesnp_TITN_gt_MAF_1_perc_finemap.log;1645

chmod 777 *
dos2unix samplesnp_TITN_gt_MAF_1_perc.z
dos2unix samplesnp_BAG3_gt_MAF_1_perc.z

#################
## Run FINEMAP ##
#################
## >>>>>>>>>>> Assuming 1 causal SNP
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/
# Results are moved here: /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP_results/

./finemap_v1.4.1_x86_64 --sss --log --n-causal-snps 1 --in-files samplesnp_TITN_gt_MAF_1_perc
# - Number of GWAS samples           : 1645
# - Number of SNPs                   : 947
# - Prior-Pr(# of causal SNPs is k)  :
#   (0 -> 0)
#    1 -> 1
# - 947 configurations evaluated
# - Computing causal SNP statistics  : done!
# - Regional SNP heritability        : 0.00173 (SD: 0.00161 ; 95% CI: [2.28e-05,0.00595])
# - Log10-BF of >= one causal SNP    : 0.0354
# - Post-expected # of causal SNPs   : 1
# - Post-Pr(# of causal SNPs is k)   :
#   (0 -> 0)
#    1 -> 1

./finemap_v1.4.1_x86_64 --sss --log --n-causal-snps 1 --in-files samplesnp_BAG3_gt_MAF_1_perc
# - Number of GWAS samples           : 1645
# - Number of SNPs                   : 1594
# - Prior-Pr(# of causal SNPs is k)  :
#   (0 -> 0)
#    1 -> 1
# - 1594 configurations evaluated
# - Computing causal SNP statistics  : done!
# - Regional SNP heritability        : 0.00209 (SD: 0.00171 ; 95% CI: [0.000105,0.00645])
# - Log10-BF of >= one causal SNP    : -0.00704
# - Post-expected # of causal SNPs   : 1
# - Post-Pr(# of causal SNPs is k)   :
#   (0 -> 0)
#    1 -> 1

## >>>>>>>>>>> Configuring each SNP of interest
# Note: log10bf column contains the log10 Bayes factors. The Bayes factor quantifies the evidence that the i-{th} SNP is causal with log10 Bayes factors greater than 2 reporting considerable evidence
./finemap_v1.4.1_x86_64 --config --in-files samplesnp_TITN_gt_MAF_1_perc --dataset 1 --rsids chr2:178562809:T:C,chr10:119670121:T:C
# - SNPs in causal configuration     : chr2:178562809:T:C
# ------------------------------------
# - Log-Likelihood of configuration  : -2330.787730
# - ML causal effect estimates       : -0.101077
# - SE of ML causal effect estimates : 0.0421795
# - P-values for causal effects      : 0.0165592
# ------------------------------------
# - Log10-BF of configuration        : 0.645256
# - Posterior mean of causal effects : -0.0813063
# - Posterior SD of causal effects   : 0.0378961
# ------------------------------------
# - Regional SNP heritability        : 0.00274 (SD: 0.00217 ; 95% CI: [8.41e-05,0.00823])



cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test
gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat  --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.00017925 --cojo-wind 100000 --out TITN_cojo
gcta64  --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat  --chr 10 --maf 0.01 --cojo-file samplesnp_BAG3_gt_MAF_1_perc_vats.ma --cojo-slct --cojo-p 0.000103976 --cojo-wind 100000 --out BAG3_cojo



module load gcta
gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-p 0.00017925 --cojo-cond cond.TTN.snplist --out cond.TTN.snplist.out
gcta64  --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat --chr 10 --maf 0.01 --cojo-file samplesnp_BAG3_gt_MAF_1_perc_vats.ma --cojo-p 0.000103976 --cojo-cond cond.BAG3.snplist --out cond.BAG3.snplist.out

# Results saved in /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test/cond.BAG3_TTN_snplist.out.cma



## Email: 08/12/2022  "I think you used summary results from SJLIFE only. Can you please re-run the analysis with summary statistics from meta-analysis?"
# USE rscript: TITN_BAG3_FINEMAP_v3_for_meta_summary_stat_email_08_12_2022.R to create samplesnp_TITN_gt_MAF_1_perc_vars_meta.ma and samplesnp_BAG3_gt_MAF_1_perc_vats_meta.ma files

# <cond.TTN.snplist>
# chr2:178562809:T:C

# <cond.BAG3.snplist>
# chr10:119670121:T:C

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test
module load gcta
gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat  --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars_meta.ma --cojo-slct --cojo-p 0.00017925 --cojo-wind 100000 --out TITN_cojo_meta
gcta64  --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat  --chr 10 --maf 0.01 --cojo-file samplesnp_BAG3_gt_MAF_1_perc_vats_meta.ma --cojo-slct --cojo-p 0.000103976 --cojo-wind 100000 --out BAG3_cojo_meta




gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars_meta.ma --cojo-p 0.00017925 --cojo-cond cond.TTN.snplist --out cond.TTN.snplist.out_meta_stat
gcta64  --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat --chr 10 --maf 0.01 --cojo-file samplesnp_BAG3_gt_MAF_1_perc_vats_meta.ma --cojo-p 0.000103976 --cojo-cond cond.BAG3.snplist --out cond.BAG3.snplist.out_meta_stat


# Results saved in /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test/cond.BAG3_TTN_snplist.out.cma_meta.xlxs 



########################################################
## Extract lines from Annotated VCF for TITN and BAG3 ##
########################################################
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/
# head -1 ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > TITN_BAG3_region_all/chr2_header
# awk '$2 ~ /^178/' ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > TITN_BAG3_region_all/chr2_178_list.txt
# cat chr2_header chr2_178_list.txt > chr2_178_list_ALL.txt

# head -1 ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > TITN_BAG3_region_all/chr10_header
# awk '$2 ~ /^119/' ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > TITN_BAG3_region_all/chr10_119_list.txt
# cat chr10_header chr10_119_list.txt > chr10_119_list_ALL.txt


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
head -1 chr2.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > ./TTN_BAG3//chr2_header_08_27_2024
awk '$2 >= 178312809 && $2 <= 178812809' chr2.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > ./TTN_BAG3/chr2_178_list_08_27_2024.txt
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3
cat chr2_header_08_27_2024 chr2_178_list_08_27_2024.txt > chr2_178_list_ALL_08_27_2024.txt

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
head -1 chr10.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > ./TTN_BAG3/chr10_header_08_27_2024
awk '$2 >= 119420121 && $2 <= 119920121'  chr10.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > ./TTN_BAG3/chr10_119_list_08_27_2024.txt
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3
cat chr10_header_08_27_2024 chr10_119_list_08_27_2024.txt > chr10_119_list_ALL_08_27_2024.txt

# The extract VCF using extract.sh with bcftools
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3

for files in *.vcf; do
plink --vcf ${files} --double-id --vcf-half-call m --keep-allele-order --threads 2 --make-bed --out plink_out/${files}.sjlife.max.maf.0.01.dat
done


cat *.bim| sort -V| uniq > all_var.sjlife.0.01.maf.txt

# check for common variants in CCSS_exp with 


# Extract overlapping P/LP
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations

plink --bfile  ../ccss_exp --bim ../ccss_exp_edited.bim --extract ../TITN.1.per.maf.gnomad.ALL.and.NFE_vars.list --max-maf 0.01 --make-bed --out TITN_ccss_exp.0.01maf
plink --bfile  ../ccss_exp --bim ../ccss_exp_edited.bim --extract ../BAG3.1.per.maf.gnomad.ALL.and.NFE_vars.list --max-maf 0.01 --make-bed --out BAG3_ccss_exp.0.01maf


plink --bfile  ../sjlife --extract ../TITN.1.per.maf.gnomad.ALL.and.NFE_vars.list --max-maf 0.01 --make-bed --out TITN_sjlife.0.01maf
plink --bfile  ../sjlife --extract ../BAG3.1.per.maf.gnomad.ALL.and.NFE_vars.list --max-maf 0.01 --make-bed --out BAG3_sjlife.0.01maf

# After extracting the common variants with maf < 0.01 in sjlife and ccss that are also rare in gnomad_ALL and gnomad_NFE, I am now extracting the overlapping variants (overlapping_vars_in_sjlife_ccss_exp_maf_0.01)

plink --bfile  TITN_ccss_exp.0.01maf --extract TITN_overlapping_vars_in_sjlife_ccss_exp_maf_0.01 --make-bed --out TITN_ccss_exp.0.01maf_overlapping
plink --bfile  BAG3_ccss_exp.0.01maf --extract BAG3_overlapping_vars_in_sjlife_ccss_exp_maf_0.01 --make-bed --out BAG3_ccss_exp.0.01maf_overlapping


plink --bfile  TITN_sjlife.0.01maf --extract TITN_overlapping_vars_in_sjlife_ccss_exp_maf_0.01 --make-bed --out TITN_sjlife.0.01maf_overlapping
plink --bfile  BAG3_sjlife.0.01maf --extract BAG3_overlapping_vars_in_sjlife_ccss_exp_maf_0.01 --make-bed --out BAG3_sjlife.0.01maf_overlapping

################################
## Now extract variants for:  ##
## 1. Clinvar and LoF         ##
## 2. Clinvar, LoF and REVEL  ##
################################

ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt .

# CCSS_exp
plink --bfile TITN_ccss_exp.0.01maf_overlapping --extract TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt --make-bed --out TITN_ccss_exp.0.01maf_overlapping.clinvar.lof
plink --bfile TITN_ccss_exp.0.01maf_overlapping --extract TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt --make-bed --out TITN_ccss_exp.0.01maf_overlapping.clinvar.lof.REVEL
plink --bfile BAG3_ccss_exp.0.01maf_overlapping --extract BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt --make-bed --out BAG3_ccss_exp.0.01maf_overlapping.clinvar.lof
plink --bfile BAG3_ccss_exp.0.01maf_overlapping --extract BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt --make-bed --out BAG3_ccss_exp.0.01maf_overlapping.clinvar.lof.REVEL

# Sjlife
plink --bfile TITN_sjlife.0.01maf_overlapping --extract TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt --make-bed --out TITN_sjlife.0.01maf_overlapping.clinvar.lof
plink --bfile TITN_sjlife.0.01maf_overlapping --extract TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt --make-bed --out TITN_sjlife.0.01maf_overlapping.clinvar.lof.REVEL
plink --bfile BAG3_sjlife.0.01maf_overlapping --extract BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt --make-bed --out BAG3_sjlife.0.01maf_overlapping.clinvar.lof
plink --bfile BAG3_sjlife.0.01maf_overlapping --extract BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt --make-bed --out BAG3_sjlife.0.01maf_overlapping.clinvar.lof.REVEL

##############
## Recode A ##
##############
# CCSS_exp

plink --bfile TITN_ccss_exp.0.01maf_overlapping --extract TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt --recode A --out TITN_ccss_exp.0.01maf_overlapping.clinvar.lof_recodeA
plink --bfile TITN_ccss_exp.0.01maf_overlapping --extract TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt --recode A  --out TITN_ccss_exp.0.01maf_overlapping.clinvar.lof.REVEL_recodeA
plink --bfile BAG3_ccss_exp.0.01maf_overlapping --extract BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt --recode A  --out BAG3_ccss_exp.0.01maf_overlapping.clinvar.lof_recodeA
plink --bfile BAG3_ccss_exp.0.01maf_overlapping --extract BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt --recode A  --out BAG3_ccss_exp.0.01maf_overlapping.clinvar.lof.REVEL_recodeA

# Sjlife
plink --bfile TITN_sjlife.0.01maf_overlapping --extract TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt --recode A --out TITN_sjlife.0.01maf_overlapping.clinvar.lof_recodeA
plink --bfile TITN_sjlife.0.01maf_overlapping --extract TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt --recode A --out TITN_sjlife.0.01maf_overlapping.clinvar.lof.REVEL_recodeA
plink --bfile BAG3_sjlife.0.01maf_overlapping --extract BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt --recode A --out BAG3_sjlife.0.01maf_overlapping.clinvar.lof_recodeA
plink --bfile BAG3_sjlife.0.01maf_overlapping --extract BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt --recode A --out BAG3_sjlife.0.01maf_overlapping.clinvar.lof.REVEL_recodeA



############################################################################################################################
## Repeat analysis for CCSS expansion (with CCM) excluding samples within 5 years of elapsed years from cancers diagnosis ##
############################################################################################################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3

# samples to remove Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_exp.young.CA.to.remove.txt

# Run GWAS: run_gwas_ccss_with_elapsed_age_5yrs.sh
# Run meta-analysis: 




## On 02/20/2023, Yadav asked me to concat all three datasets (sjlife, CCSS_org and CCSS_exp) and run the assocation analysis on the merged data; see Concatenated analysis on 02/20/2023 for phenotype file in 3.demographics_v2.R
plink --bfile sjlife --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out sjlife_to_concat
plink --bfile ccss_org --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out ccss_org_to_concat
plink --bfile ccss_exp --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out ccss_exp_to_concat

## Merge the above three
plink --bfile sjlife_to_concat --merge-list merge_list.txt --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples
<merged_list.txt>
ccss_org_to_concat
ccss_exp_to_concat

## Calculate frequency
plink --bfile sjlife_ccss_org_ccss_exp_samples --freq --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples_freq_out



# Harmonize alleles using merge association files
plink --bfile sjlife_to_concat --update-name harmonize.txt --make-bed --out sjlife_to_concat_updated
plink --bfile ccss_org_to_concat --update-name harmonize.txt --make-bed --out ccss_org_to_concat_updated
plink --bfile ccss_exp_to_concat --update-name harmonize.txt --make-bed --out ccss_exp_to_concat_updated

## Merge the above three
plink --bfile sjlife_to_concat_updated --merge-list merge_list.txt --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples_updated
<merged_list.txt>
ccss_org_to_concat_updated
ccss_exp_to_concat_updated

## Calculate frequency
plink --bfile sjlife_ccss_org_ccss_exp_samples_updated --freq --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples_updated_freq_out


## Run association analysis
 plink \
 --allow-no-sex \
 --bfile sjlife_ccss_org_ccss_exp_samples_updated \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,cohort_two \
 --out sjlife_ccss_org_ccss_exp_samples_results




## Run association analysis with ccss_org and ccss_exp
plink --bfile ccss_org_to_concat_updated --bmerge ccss_exp_to_concat_updated --keep-allele-order --make-bed --out merged_ccss


 # plink \
 # --allow-no-sex \
 # --bfile merged_ccss \
 # --maf 0.01 \
 # --hwe 1e-06 \
 # --logistic hide-covar \
 # --snp assoc_test_two_vars.txt \
 # --ci 0.95 \
 # --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 # --pheno-name CMP \
 # --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 # --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 # --out ccss_org_ccss_exp_vars


plink \
 --allow-no-sex \
 --bfile merged_ccss \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out ccss_org_ccss_exp


plink \
 --allow-no-sex \
 --bfile ccss_org_to_concat_updated \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out ccss_org

plink \
 --allow-no-sex \
 --bfile ccss_exp_to_concat_updated \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out ccss_exp


plink \
 --allow-no-sex \
 --bfile sjlife_to_concat_updated \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out sjlife_results


## AFR
# plink 
#   --bfile sjlife
#   --extract common_missense_variants.txt
#   --keep pheno/sjlife_ttn_bag3_afr.pheno
#   --make-bed
#   --out sjlife_afr_ttn_bag3

tail -n +2 common_missense_variants.txt | cut -f1> extract_common_missense_vars
tail -n +2 pheno/sjlife_ttn_bag3_afr.pheno | cut -f1,2 > afr_samples

plink --bfile sjlife --update-name harmonize.txt --make-bed --out sjlife_afr_to_concat_updated
plink --bfile sjlife_afr_to_concat_updated --extract extract_common_missense_vars --keep afr_samples --make-bed --out sjlife_afr_ttn_bag3_AN

plink --bfile sjlife_afr_ttn_bag3_AN

plink \
  --1 \
  --allow-no-sex \
  --bfile sjlife_afr_ttn_bag3_AN \
  --ci 0.95 \
  --covar pheno/sjlife_ttn_bag3_afr.pheno \
  --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --hwe 1e-06 \
  --logistic hide-covar \
  --out sjlife_results_afr_AN \
  --pheno pheno/sjlife_ttn_bag3_afr.pheno \
  --pheno-name CMP

plink --bfile sjlife_afr_ttn_bag3_AN --freq --keep afr_samples --out sjlife_results_afr_freq_AN










## LD check (based on Kateryna's email on 03/17/2023)
# 2 (window size kb, needs to be more than 2) 1 (step size, it could be 1 or 2 in such a small dataset) and 0.8 (it is r2 threshold), generated  plink.prune.in will contain variants that pass the threshold.
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3
plink --bfile ttn_significant_ld --indep-pairwise 2 1 0.8


# ## Keep PED with plink.prune.in variants
# plink --bfile ttn_significant_ld --extract plink.prune.in --keep-allele-order --recode --out haplotype_input
plink --bfile ttn_significant_ld --extract plink.prune.in --keep-allele-order --recodeA --out haplotype_input

# plink --bfile ttn_significant_ld --recode --out haplotype_input
# plink --file haplotype_input --hap plink.prune.in 

## In plink 1.70
# plink --bfile ttn_significant_ld --hap plink.prune.in --noweb


# ## Since got this error in plink 1.90: Error: The --hap... family of flags has not been reimplemented in PLINK 1.9 due; I am using beagle
#  plink --bfile ttn_significant_ld --recode beagle --out ttn_significant_ld_beagle

## Run convert_raw_to_phasing_v2.R script, then run PHASE
 PHASE haplotype_input_edited.txt haplotype_phase.out


#  ## Repeat with r2 0.2
#  ## LD check (based on Kateryna's email on 03/17/2023)
# # 2 (window size kb, needs to be more than 2) 1 (step size, it could be 1 or 2 in such a small dataset) and 0.2 (it is r2 threshold), generated  plink.prune.in will contain variants that pass the threshold.
# # also, just to note, this is a random pruning, so I guess we will need to make sure that after pruning our snp of interest remains in the list (we could just repeat the pruning with different step if it is removed - instead of 1 use 2 or 3 for example)
# plink --bfile ttn_significant_ld --indep-pairwise 2 1 0.2 --out prune_0.02

# # ## Keep PED with plink.prune.in variants
# # plink --bfile ttn_significant_ld --extract plink.prune.in --keep-allele-order --recode --out haplotype_input
# plink --bfile ttn_significant_ld --extract prune_0.02.prune.in --keep-allele-order --recodeA --out haplotype_input_0.2
# PHASE -c haplotype_input_edited_0.2.txt haplotype_phase_0.2.out



