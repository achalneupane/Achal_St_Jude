#############################
## Common variant analysis ##
#############################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3
## On 02/20/2023, Combined all three datasets (sjlife, CCSS_org and CCSS_exp) and run the assocation analysis on the merged data; see Concatenated analysis on 02/20/2023 for phenotype file in 3.demographics_v2.R
plink --bfile sjlife --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out sjlife_to_concat
plink --bfile ccss_org --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out ccss_org_to_concat
plink --bfile ccss_exp --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out ccss_exp_to_concat

## Merge the above three
plink --bfile sjlife_to_concat --merge-list merge_list.txt --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples

## Calculate frequency
plink --bfile sjlife_ccss_org_ccss_exp_samples --freq --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples_freq_out

## Run association analysis (Combined: sjlife, ccss_exp, ccss_org)
 plink \
 --allow-no-sex \
 --bfile sjlife_ccss_org_ccss_exp_samples \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out sjlife_ccss_org_ccss_exp_samples_results


# Association analysis on individual datasets
## Run GWAS for each cohort
module load plink/1.90b
# SJLIFE
 plink \
 --allow-no-sex \
 --bfile sjlife_to_concat \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno  \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno  \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out sjlife_results_AN
# CCSS org
 plink \
 --allow-no-sex \
 --bfile ccss_org_to_concat \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno  \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out ccss_org_results_AN
# CCSS exp
 plink \
 --allow-no-sex \
 --bfile ccss_exp_to_concat \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno  \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno  \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out ccss_exp_results_AN






# [7/13 4:09 PM] Petrykey, Kateryna

# Hey Achal, I was working on the TTN manuscript., this version that you sent me is definitely much better, but we still need to work a lot on Discussion. There are a few things we need to understand in order to explain why we haven't seen the effect of rare variants reported in cancer patients by previous study. First, we need to check which variants exactly have been tested by Garcia-Pavia, P., et al. Genetic Variants Associated With Cancer Therapy-Induced Cardiomyopathy. Circulation 140, 31-41 (2019). if we have tested the same variants, and if not, we will need to investigate why these variants were not available in our dataset. Additionally, we will need to look for a specific location of the variants included in our analysis. It has been consistently demonstrated that truncated variants linked to cardiomyopathy were located on A-band of the protein, so we should probably include a figure with the location of the variants.

# [7/13 4:14 PM] Petrykey, Kateryna

# Also, several studies reported that male carriers were more prone to cardiomyopathy than females, so we may want to check this too in our study, it will probably will not make a lot of sense for rare variants analysis, but we could definitely check for common and echo parameters. Please let me know what you think, thanks a lot.

# In this command, I've added the --interaction gender option, which will include an interaction term between gender and the carrier status (presumably specified in your phenotype file under "CMP"). This interaction term will allow you to assess whether the effect of being a carrier on cardiomyopathy risk differs between males and females.

# After running this modified command, you can look at the results in the output file (sjlife_ccss_org_ccss_exp_samples_results_gender_interaction.log) to determine whether male carriers have a different risk of cardiomyopathy compared to female carriers. Look for the interaction term's p-value to assess the statistical significance of this difference. If the p-value is significant, it suggests that there is a differential effect of being a carrier on cardiomyopathy risk between males and females.

## 0 is male; 1 is female

# awk '$6 == 1 || NR == 1' pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno > pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_female.pheno
# awk '$6 == 0 || NR == 1' pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno > pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_male.pheno



# ## Male
#  plink \
#  --allow-no-sex \
#  --bfile sjlife_ccss_org_ccss_exp_samples \
#  --maf 0.01 \
#  --hwe 1e-06 \
#  --logistic hide-covar \
#  --ci 0.95 \
#  --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_male.pheno \
#  --pheno-name CMP \
#  --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_male.pheno \
#  --covar-name agedx,agelstcontact,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
#  --out sjlife_ccss_org_ccss_exp_samples_results_male


# ## Female
#  plink \
#  --allow-no-sex \
#  --bfile sjlife_ccss_org_ccss_exp_samples \
#  --maf 0.01 \
#  --hwe 1e-06 \
#  --logistic hide-covar \
#  --ci 0.95 \
#  --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_female.pheno \
#  --pheno-name CMP \
#  --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_female.pheno \
#  --covar-name agedx,agelstcontact,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
#  --out sjlife_ccss_org_ccss_exp_samples_results_female

## Run association analysis with ccss_org and ccss_exp
plink --bfile ccss_org_to_concat --bmerge ccss_exp_to_concat --make-bed --out merged_ccss


 plink \
 --allow-no-sex \
 --bfile merged_ccss \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --extract assoc_test_two_vars.txt \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out sccss_org_ccss_exp_vars




 ## Analysis of chr2_178562809_T_C and chr10:119670121:T:C carriers:
 # Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\combined_cohorts\cardiomyopathy_gwas_echo\echo_dosage_and_gender_analysis\TTN_BAG3_carrier_heartRT_anthra\top_TTN_BAG3_SNPs_wrt_carrier_status_in_sjlife_ccss_exp_ccss_org.R
plink --bfile sjlife_ccss_org_ccss_exp_samples --out sjlife_ccss_org_ccss_exp_chr2_178562809_T_C --recode A --snp chr2:178562809:T:C
plink --bfile sjlife_ccss_org_ccss_exp_samples --out sjlife_ccss_org_ccss_exp_chr10_119670121_T_C --recode A --snp chr10:119670121:T:C
