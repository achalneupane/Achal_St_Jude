#############################
## Common variant analysis ##
#############################
## On 02/20/2023, Combined all three datasets (sjlife, CCSS_org and CCSS_exp) and run the assocation analysis on the merged data; see Concatenated analysis on 02/20/2023 for phenotype file in 3.demographics_v2.R
plink --bfile sjlife --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out sjlife_to_concat
plink --bfile ccss_org --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out ccss_org_to_concat
plink --bfile ccss_exp --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out ccss_exp_to_concat

## Merge the above three
plink --bfile sjlife_to_concat --merge-list merge_list.txt --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples

## Calculate frequency
plink --bfile sjlife_ccss_org_ccss_exp_samples --freq --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples_freq_out

## Run association analysis
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