# # Additionally, .z file shoul have maf in specific form. I was getting this error: "Error : Expected a floating point value greater than 0.0 and smaller or equal than 0.5 in line 3 column 6 of file 'samplesnp.z'!"

## Extract TITN and BAG3 variants with MAF > 1%; separately.
plink --bfile samplesnp.dat --extract samplesnp_TITN_gt_MAF_1_perc_vars.list --keep-allele-order --make-bed --out samplesnp_TITN_gt_MAF_1_perc_vars.dat

# RUN GEC analysis to calculates the number of independent variants in datasets prepared above <http://pmglab.top/gec/#/download>; wget http://pmglab.top/gec/data/archive/v0.2/gecV0.2.zip
module load java
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/gec
ln -s ../samplesnp_TITN_gt_MAF_1_perc_vars.dat* .


java -jar -Xmx1g gec.jar --var-gec --effect-number --plink-binary samplesnp_TITN_gt_MAF_1_perc_vars.dat --maf 0.01  --genome --out samplesnp_TITN_gt_MAF_1_perc_vars.dat_GEC_test1

# calculate Bonferroni corrected P-value
(base) [aneupane@splprhpc09 gec]$ cat samplesnp_TITN_gt_MAF_1_perc_vars.dat_GEC_test1.sum
Observed_Number Effective_Number        Effective_Ratio Suggestive_P_Value      Significant_P_Value     Highly_Significant_P_Value
947     278.94  0.29    3.58E-3 1.79E-4 3.58E-6

BF_P_TITN = 0.00017925 = 0.05/278.94 


## Run conditional joint analysis; https://gcta.freeforums.net/thread/178/conditional-joint-analysis-using-summary
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test
ln -s ../samplesnp_TITN_gt_MAF_1_perc_vars.dat* .


module load gcta


<R>
## First only keep those with maf >1 %; then separate the summary statistics for TITN and BAG3 and 
gwas.dat.maf.gt.1perc <- gwas.dat[!gwas.dat$maf < 0.01,]
gwas.dat.maf.gt.1perc.TITN <- gwas.dat.maf.gt.1perc[grepl("2", gwas.dat.maf.gt.1perc$chromosome),]
write.table(as.data.frame(gwas.dat.maf.gt.1perc.TITN$rsid), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp_TITN_gt_MAF_1_perc_vars_meta.list", quote = F, row.names = F, col.names = F, sep = " ")
# After extracting these variants from plink, I am re-ordering the summary stat
# to match the order of varinats in plink. This is to ensure LD and summary stat
# have same variant order
samplesnp_TITN_gt_MAF_1_perc_vars.dat.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp_TITN_gt_MAF_1_perc_vars.dat.bim")
gwas.dat.maf.gt.1perc.TITN <- gwas.dat.maf.gt.1perc.TITN[match(samplesnp_TITN_gt_MAF_1_perc_vars.dat.bim$V2, gwas.dat.maf.gt.1perc.TITN$rsid),]
table(gwas.dat.maf.gt.1perc.TITN$rsid == samplesnp_TITN_gt_MAF_1_perc_vars.dat.bim$V2)
# TRUE 
# 947 
write.table(gwas.dat.maf.gt.1perc.TITN, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp_TITN_gt_MAF_1_perc_meta.z", quote = F, row.names = F, col.names = T, sep = " ")
## Input summary stat file for COJO analysis 
TITN.COJO <- cbind.data.frame(SNP=gwas.dat.maf.gt.1perc.TITN$rsid, A1=gwas.dat.maf.gt.1perc.TITN$allele1, A2=gwas.dat.maf.gt.1perc.TITN$allele2, 
                              freq=gwas.dat.maf.gt.1perc.TITN$maf, b=gwas.dat.maf.gt.1perc.TITN$beta, se=gwas.dat.maf.gt.1perc.TITN$se, p=gwas.dat.maf.gt.1perc.TITN$P)
TITN.COJO$N <- 1645
## COJO files
write.table(TITN.COJO, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test/samplesnp_TITN_gt_MAF_1_perc_vars_meta.ma", quote = F, row.names = F, col.names = F, sep = " ")


less -S samplesnp_TITN_gt_MAF_1_perc_vars.ma
chr2:178313079:G:A A G 0.0845 0.231111720963387 0.190410329535685 0.213 1645
chr2:178313658:T:C T C 0.01459 -1.10866262452161 0.733656312073577 0.13 1645
chr2:178313672:C:A C A 0.01459 -1.10866262452161 0.733656312073577 0.13 1645
chr2:178313675:A:G G A 0.08399 0.2390169004705 0.189157485280813 0.201 1645
chr2:178313779:A:G A G 0.01459 -1.10866262452161 0.733656312073577 0.13 1645


gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat  --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.00017925 --cojo-wind 100000 --out TITN_cojo


## Test if GCTA options are valid with p 0.01 or even 0.05
# gcta64  --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat  --chr 2 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.05 --cojo-wind 100000 --out TITN_cojo


## Since gcta analysis did not find any independent SNPs, which implies we can assume 1 causal variant within each of the TTN and BAG3 locus



# Create LD file
plink --r2 --bfile samplesnp_TITN_gt_MAF_1_perc_vars.dat --matrix --out samplesnp_TITN_gt_MAF_1_perc

<samplesnp.z>
rsid chromosome position allele1 allele2 maf beta se
chr2:178313079:G:A 2 178313079 A G 0.0845 0.231111720963387 0.190410329535685
chr2:178313658:T:C 2 178313658 T C 0.01459 -1.10866262452161 0.733656312073577
chr2:178313672:C:A 2 178313672 C A 0.01459 -1.10866262452161 0.733656312073577
chr2:178313675:A:G 2 178313675 G A 0.08399 0.2390169004705 0.189157485280813



## <samplesnp_TITN_gt_MAF_1_perc>
# z;ld;snp;config;cred;log;n_samples
# samplesnp_TITN_gt_MAF_1_perc.z;samplesnp_TITN_gt_MAF_1_perc.ld;samplesnp_TITN_gt_MAF_1_perc.snp;samplesnp_TITN_gt_MAF_1_perc.config;samplesnp_TITN_gt_MAF_1_perc.cred;samplesnp_TITN_gt_MAF_1_perc_finemap.log;1645


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




