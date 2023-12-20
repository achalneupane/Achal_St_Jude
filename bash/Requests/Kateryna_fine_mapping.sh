# VCF to plink
plink --double-id --keep-allele-order --make-bed --out /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR//MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed --threads 4 --vcf MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz --vcf-half-call m

# calculate frequency
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas
# plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed --keep-allele-order --keep pheno/sjlife_eur_dox_only_pcs.pheno --freq --out MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_freq_993samples
# awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_freq_993samples.frq > MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_freq_993samples.frq_edited1
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed --keep-allele-order --freq --out MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_freq
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_freq.frq > MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_freq.frq_edited1



<R: finemapping.R>
# once you run R code, run :
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR//MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed --extract samplesnp_gt_MAF_1_perc_vars_meta.list --keep ../pheno/sjlife_eur_dox_only_pcs.pheno --keep-allele-order --make-bed --out chr16_finemap_plink
#####


# # Additionally, .z file shoul have maf in specific form. I was getting this error: "Error : Expected a floating point value greater than 0.0 and smaller or equal than 0.5 in line 3 column 6 of file 'samplesnp.z'!"




# RUN GEC analysis to calculates the number of independent variants in datasets prepared above <http://pmglab.top/gec/#/download>; wget http://pmglab.top/gec/data/archive/v0.2/gecV0.2.zip
module load java
java -jar -Xmx1g /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/gec/gec.jar --var-gec --effect-number --plink-binary chr16_finemap_plink --maf 0.01  --genome --out chr16_finemap_plink_GEC_test1

# calculate Bonferroni corrected P-value
(base) [aneupane@splprhpc09 gec]$ cat chr16_finemap_plink_GEC_test1.sum
Observed_Number Effective_Number        Effective_Ratio Suggestive_P_Value      Significant_P_Value     Highly_Significant_P_Value
847     172.42  0.2     5.80E-3 2.90E-4 5.80E-6

BF_P_TITN = 0.00028999 = 0.05/172.42


## Run conditional joint analysis; https://gcta.freeforums.net/thread/178/conditional-joint-analysis-using-summary
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test
ln -s ../chr16_finemap_plink* .


module load gcta





less -S samplesnp_gt_MAF_1_perc_vars.ma
chr2:178313079:G:A A G 0.0845 0.231111720963387 0.190410329535685 0.213 1645
chr2:178313658:T:C T C 0.01459 -1.10866262452161 0.733656312073577 0.13 1645
chr2:178313672:C:A C A 0.01459 -1.10866262452161 0.733656312073577 0.13 1645
chr2:178313675:A:G G A 0.08399 0.2390169004705 0.189157485280813 0.201 1645
chr2:178313779:A:G A G 0.01459 -1.10866262452161 0.733656312073577 0.13 1645

# the command is using GCTA to perform a COJO analysis on chromosome 16, considering variants with a MAF above 0.01. The analysis includes selecting independent significant SNPs based on a specified p-value threshold within a defined genomic window. The results will be saved with the prefix "cojo."
gcta64  --bfile chr16_finemap_plink  --chr 16 --maf 0.01 --cojo-file samplesnp_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.00028999 --cojo-wind 100000 --out cojo

# Create LD file
plink --r2 --bfile chr16_finemap_plink --matrix --out samplesnp_gt_MAF_1_perc
## test
plink --r2 --bfile chr16_finemap_plink --out samplesnp_gt_MAF_1_perc_test_ld

(base) [aneupane@noderome105 finemap_kateryna]$ head samplesnp_gt_MAF_1_perc.z
rsid chromosome position allele1 allele2 maf beta se P
chr16:25362303:A:G 16 25362303 G A 0.12 0.0363319292473902 0.1723 0.831
chr16:25362595:A:G 16 25362595 G A 0.169 -0.270890509811944 0.1634 0.09744
chr16:25363139:C:A 16 25363139 A C 0.4505 0.180653499693258 0.1152 0.1175
chr16:25364191:G:T 16 25364191 T G 0.4507 0.179818426575836 0.1153 0.1182
chr16:25364505:A:G 16 25364505 G A 0.1469 0.0779615414697118 0.1534 0.611
chr16:25364735:A:G 16 25364735 G A 0.08694 -0.0769890415281357 0.2651 0.7714




## <samplesnp_gt_MAF_1_perc>
# z;ld;snp;config;cred;log;n_samples
# samplesnp_gt_MAF_1_perc.z;samplesnp_gt_MAF_1_perc.ld;samplesnp_gt_MAF_1_perc.snp;samplesnp_gt_MAF_1_perc.config;samplesnp_gt_MAF_1_perc.cred;samplesnp_gt_MAF_1_perc_finemap.log;1645


# <samplesnp_gt_MAF_1_perc>
z;ld;snp;config;cred;log;n_samples
samplesnp_gt_MAF_1_perc.z;samplesnp_gt_MAF_1_perc.ld;samplesnp_gt_MAF_1_perc.snp;samplesnp_gt_MAF_1_perc.config;samplesnp_gt_MAF_1_perc.cred;samplesnp_gt_MAF_1_perc_finemap.log;993




#################
## Run FINEMAP ##
#################
## >>>>>>>>>>> Assuming 1 causal SNP
/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64
# Results are moved here: /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP_results/

/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 --sss --log --in-files samplesnp_gt_MAF_1_perc
# ------------
# FINE-MAPPING (1/1)
# ------------
# - GWAS summary stats               : samplesnp_gt_MAF_1_perc.z
# - SNP correlations                 : samplesnp_gt_MAF_1_perc.ld
# - Causal SNP stats                 : samplesnp_gt_MAF_1_perc.snp
# - Causal configurations            : samplesnp_gt_MAF_1_perc.config
# - Credible sets                    : samplesnp_gt_MAF_1_perc.cred
# - Log file                         : samplesnp_gt_MAF_1_perc_finemap.log_sss
# - Reading input                    : done!

# - Updated prior SD of effect sizes : 0.05 0.0644 0.083 0.107

# - Number of GWAS samples           : 993
# - Number of SNPs                   : 859
# - Prior-Pr(# of causal SNPs is k)  :
#   (0 -> 0)
#    1 -> 0.583
#    2 -> 0.291
#    3 -> 0.097
#    4 -> 0.0242
#    5 -> 0.00482
# - 7830 configurations evaluated (0.369/100%) : converged after 369 iterations
# - Computing causal SNP statistics  : done!
# - Regional SNP heritability        : 0.0255 (SD: 0.00914 ; 95% CI: [0.0103,0.0457])
# - Log10-BF of >= one causal SNP    : 3.47
# - Post-expected # of causal SNPs   : 1.41
# - Post-Pr(# of causal SNPs is k)   :
#   (0 -> 0)
#    1 -> 0.608
#    2 -> 0.376
#    3 -> 0.0151
#    4 -> 0
#    5 -> 0
# - Writing output                   : done!
# - Run time                         : 0 hours, 0 minutes, 1 seconds


## >>>>>>>>>>> Configuring each SNP of interest
# Note: log10bf column contains the log10 Bayes factors. The Bayes factor quantifies the evidence that the i-{th} SNP is causal with log10 Bayes factors greater than 2 reporting considerable evidence
# Evaluate a single causal configuration without performing shotgun stochastic search (--config)

# # https://www.biorxiv.org/content/biorxiv/early/2022/10/24/2022.10.21.513123.full.pdf
# First, these methods determine
# a posterior inclusion probability (PIP) for each variant, quantifying the probability that the variant
# is causal under the model, which can reflect uncertainty due to LD. For example, two variants in
# perfect LD and harboring a strong association with the phenotype may each have PIP 50%,
# representing confidence that one but not likely both variants are causal. Second, these methods
# incorporate assumptions about genetic architecture -- namely, the relative probabilities of
# different numbers of and configurations of causal SNPs, as reflected by a Bayesian prior -- to
# improve statistical power for identifying high-confidence variants.
# Bayesian fine-mapping methods are correctly calibrated when the PIPs accurately reflect the
# true proportions of causal variants, e.g. 9 out of 10 variants having PIP 90% are truly causal for
# the trait. Calibration is ensured when the linear model for genetic effects and Bayesian prior for
# genetic architecture across loci are both correctly specified, and accurate calibration has also
# been demonstrated empirically in simulations to be robust under mild model misspecifications2

#############
## Results ##
#############

/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP_results/without_assumption_4_causal_SNP

base) [aneupane@splprhpc09 without_assumption_4_causal_SNP]$ head  samplesnp_TITN_gt_MAF_1_perc.config
rank config prob log10bf odds k prob_norm_k h2 h2_0.95CI mean sd
1 chr2:178395655:A:G 0.014982 1.18733 1 1 0.014982 0.00396585 0.000306971,0.0110977 -0.084899 0.0317607
2 chr2:178469905:C:T 0.00803431 0.916714 1.86475 1 0.00803431 0.00335447 0.000174501,0.00935737 0.33325 0.137379
3 chr2:178763642:G:T 0.00729498 0.87479 2.05373 1 0.00729498 0.0032598 0.00012164,0.0087938 -0.101792 0.0426758
4 chr2:178617620:CTAAG:C 0.00666275 0.835419 2.24862 1 0.00666275 0.0031709 0.000143807,0.00944238 0.0902174 0.0384472
5 chr2:178519690:T:C 0.00645174 0.821443 2.32216 1 0.00645174 0.00313934 5.44437e-05,0.00889824 -0.0894485 0.0383466
6 chr2:178516596:T:G 0.00645174 0.821443 2.32216 1 0.00645174 0.00313934 0.000112958,0.00890202 -0.0894798 0.03836
7 chr2:178599800:T:C 0.00598226 0.788631 2.5044 1 0.00598226 0.00306526 7.05721e-05,0.00877184 -0.0867095 0.0377047
8 chr2:178677480:T:A 0.00544996 0.748159 2.749 1 0.00544996 0.0029739 6.21106e-05,0.00872368 0.317734 0.140688
9 chr2:178548239:G:A 0.00544567 0.747817 2.75117 1 0.00544567 0.00297313 4.49943e-05,0.00894418 -0.217308 0.0962361


Let's interpret the columns in this file:
rank: The rank of the configuration. It indicates the order in which the configurations are ranked based on their posterior probability.
config: The configuration of the causal variant. It includes information about the chromosome, position, and alleles.
prob: The posterior probability of the configuration being the true causal variant.
log10bf: The base-10 logarithm of the Bayes factor. It quantifies the strength of evidence for each configuration.
odds: The odds ratio for the configuration being causal.
k: The number of causal variants in the configuration.
prob_norm_k: The normalized probability of the configuration, considering the number of causal variants.
h2: The estimated heritability explained by the configuration.
h2_0.95CI: The 95% confidence interval for the heritability.
mean: The mean effect size of the causal variant.
sd: The standard deviation of the effect size.

Here's an interpretation of the first few lines:
The configuration with rank 1 has a variant on chromosome 2 at position 178395655, where the reference allele is A and the alternative allele is G. The posterior probability of this configuration being the true causal variant is 0.014982. The Bayes factor is 1.18733, and the odds ratio is 1. The configuration includes 1 causal variant. The normalized probability considering the number of causal variants is 0.014982. The estimated heritability explained by this configuration is -0.084899, with a 95% confidence interval of 0.000306971 to 0.0110977. The mean effect size is -0.084899, and the standard deviation is 0.0317607.
Similarly, each subsequent line provides information about different configurations ranked by their posterior probability.
These results are indicative of FineMap's attempt to identify the most likely causal variants in the specified genomic region and provide estimates of their effects and heritability. Keep in mind that the interpretation may vary based on the specific assumptions and parameters used in the FineMap analysis.

A moderate threshold (e.g., 0.5 to 0.75) can be a reasonable compromise between sensitivity and specificity. This approach aims to balance the risk of false positives and false negatives.