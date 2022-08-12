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
mv *.cred1 *.bf1 *.config *.snp *.log_sss /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP_results/assuming_1_causal_SNP/

## >>>>>>>>>>> Without any asumptions for the number of causal SNP
./finemap_v1.4.1_x86_64 --sss --log --in-files samplesnp_TITN_gt_MAF_1_perc --dataset 1
# - Number of GWAS samples           : 1645
# - Number of SNPs                   : 947
# - Prior-Pr(# of causal SNPs is k)  :
#   (0 -> 0)
#    1 -> 0.583
#    2 -> 0.291
#    3 -> 0.097
#    4 -> 0.0242
#    5 -> 0.00482
# - 947 configurations evaluated (0.100/100%) : converged after 100 iterations
# - Computing causal SNP statistics  : done!
# - Regional SNP heritability        : 0.00173 (SD: 0.00161 ; 95% CI: [2.26e-05,0.00593])
# - Log10-BF of >= one causal SNP    : -0.199
# - Post-expected # of causal SNPs   : 1
# - Post-Pr(# of causal SNPs is k)   :
#   (0 -> 0)
#    1 -> 1
#    2 -> 0
#    3 -> 0
#    4 -> 0
#    5 -> 0

./finemap_v1.4.1_x86_64 --sss --log --in-files samplesnp_BAG3_gt_MAF_1_perc --dataset 1
# - Number of GWAS samples           : 1645
# - Number of SNPs                   : 1594
# - Prior-Pr(# of causal SNPs is k)  :
#   (0 -> 0)
#    1 -> 0.583
#    2 -> 0.291
#    3 -> 0.097
#    4 -> 0.0242
#    5 -> 0.00484
# - 1604 configurations evaluated (0.133/100%) : converged after 133 iterations
# - Computing causal SNP statistics  : done!
# - Regional SNP heritability        : 0.00212 (SD: 0.00172 ; 95% CI: [0.00013,0.00654])
# - Log10-BF of >= one causal SNP    : -0.24
# - Post-expected # of causal SNPs   : 1
# - Post-Pr(# of causal SNPs is k)   :
#   (0 -> 0)
#    1 -> 0.995
#    2 -> 0.00489
#    3 -> 0
#    4 -> 0
#    5 -> 0

mv *.cred1 *.bf1 *.config *.snp *.log_sss /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP_results/without_assumption_4_causal_SNP


## >>>>>>>>>>> Configuring each SNP of interest
# Note: log10bf column contains the log10 Bayes factors. The Bayes factor quantifies the evidence that the i-{th} SNP is causal with log10 Bayes factors greater than 2 reporting considerable evidence
./finemap_v1.4.1_x86_64 --config --in-files samplesnp_TITN_gt_MAF_1_perc --dataset 1 --rsids chr2:178562809:T:C
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


./finemap_v1.4.1_x86_64 --config --in-files samplesnp_BAG3_gt_MAF_1_perc --dataset 1 --rsids chr10:119670121:T:C
# - SNPs in causal configuration     : chr10:119670121:T:C
# ------------------------------------
# - Log-Likelihood of configuration  : -2331.429058
# - ML causal effect estimates       : -0.0896417
# - SE of ML causal effect estimates : 0.0424663
# - P-values for causal effects      : 0.0347816
# ------------------------------------
# - Log10-BF of configuration        : 0.421902
# - Posterior mean of causal effects : -0.0721079
# - Posterior SD of causal effects   : 0.0381389
# ------------------------------------
# - Regional SNP heritability        : 0.00224 (SD: 0.00197 ; 95% CI: [1.42e-05,0.0072])

# - Run time                         : 0 hours, 0 minutes, 0 seconds



## Email: 08/11/2022  "can you please also run GCTA analysis but conditioning on the SNPs in bold for both TTN and BAG3 separately?"
gcta64  --bfile test  --chr 1 --maf 0.01 --cojo-file test.ma --cojo-cond cond.bag3.snplist --out test_chr1

# <cond.TTN.snplist>
# chr2:178562809:T:C

# <cond.BAG3.snplist>
# chr10:119670121:T:C

gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat  --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.00017925 --cojo-wind 100000 --out TITN_cojo
gcta64  --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat  --chr 10 --maf 0.01 --cojo-file samplesnp_BAG3_gt_MAF_1_perc_vats.ma --cojo-slct --cojo-p 0.000103976 --cojo-wind 100000 --out BAG3_cojo


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64
module load gcta
gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-p 0.00017925 --cojo-cond cond.TTN.snplist --out cond.TTN.snplist.out
gcta64  --bfile samplesnp_BAG3_gt_MAF_1_perc_vars.dat --chr 10 --maf 0.01 --cojo-file samplesnp_BAG3_gt_MAF_1_perc_vats.ma --cojo-p 0.000103976 --cojo-cond cond.BAG3.snplist --out cond.BAG3.snplist.out


########################################################
## Extract lines from Annotated VCF for TITN and BAG3 ##
########################################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/
head -1 ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > TITN_BAG3_region_all/chr2_header
awk '$2 ~ /^178/' ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > TITN_BAG3_region_all/chr2_178_list.txt
cat chr2_header chr2_178_list.txt > chr2_178_list_ALL.txt

head -1 ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > TITN_BAG3_region_all/chr10_header
awk '$2 ~ /^119/' ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > TITN_BAG3_region_all/chr10_119_list.txt
cat chr10_header chr10_119_list.txt > chr10_119_list_ALL.txt


# The extract VCF using extract.sh with bcftools
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3

for files in *.vcf; do
plink --vcf ${files} --double-id --vcf-half-call m --keep-allele-order --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/samples_for_maf.txt --threads 2 --make-bed --max-maf 0.01 --out plink_out/${files}.sjlife.max.maf.0.01.dat
done


cat *.bim| sort -V| uniq > all_var.sjlife.0.01.maf.txt

# check for common variants in CCSS_exp with 