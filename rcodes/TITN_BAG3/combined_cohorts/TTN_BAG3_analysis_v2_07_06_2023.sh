#############################
## Common variant analysis ##
#############################
## On 02/20/2023, Combined all three datasets (sjlife, CCSS_org and CCSS_exp) and run the assocation analysis on the merged data; see Concatenated analysis on 02/20/2023 for phenotype file in 3.demographics_v2.R
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3
plink --bfile sjlife --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out sjlife_to_concat
plink --bfile ccss_org --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out ccss_org_to_concat
plink --bfile ccss_exp --keep sjlife_ccss_org_ccss_exp_samples.txt --keep-allele-order --make-bed --out ccss_exp_to_concat



# Run the Rscript to fix problematic alleles: combine_variants_from_sjlife_ccss_exp_and_ccss_org.R

## Merge the above three
plink --bfile sjlife_to_concat --merge-list merge_list.txt --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples


## extract problem variants
grep "Warning" sjlife_ccss_org_ccss_exp_samples.log | awk -F "'" '{print $2}'|awk -F "[: ]" '{print $1":"$2}'|sort -V| uniq > warnings.txt


## exclude problematic variants from merged file
plink --bfile sjlife_ccss_org_ccss_exp_samples --exclude extract_vars.txt --make-bed --out good_dataset

# ## update names in good ones
plink --bfile good_dataset --update-name update_names_good_snps.txt --make-bed --out good_dataset_updated


## extract problematic variants
plink --bfile sjlife_ccss_org_ccss_exp_samples --extract extract_vars.txt --make-bed --out extracted_dataset



plink --bfile extracted_dataset --update-alleles update_alleles.txt --make-bed --out extracted_dataset_updated_alleles


## update names
plink --bfile extracted_dataset_updated_alleles --update-name update_names.txt --make-bed --out extracted_dataset_updated_names



## merge good datasets with updated dataset
plink --bfile good_dataset_updated --bmerge  extracted_dataset_updated_names --make-bed --out sjlife_ccss_org_ccss_exp_samples_v2


## Calculate frequency
plink --bfile sjlife_ccss_org_ccss_exp_samples_v2 --freq --keep-allele-order --out sjlife_ccss_org_ccss_exp_samples_freq_out_v2

## Run association analysis
 plink \
 --allow-no-sex \
 --bfile sjlife_ccss_org_ccss_exp_samples_v2 \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out sjlife_ccss_org_ccss_exp_samples_results_v2


Final analysis results are saved as sjlife_ccss_org_ccss_exp__ttn_bag3.assoc.final_v2_07_06_2023.xlxs



## Run association analysis with ccss_org and ccss_exp
plink --bfile ccss_org_to_concat --bmerge ccss_exp_to_concat --make-bed --out merged_ccss


# Then Run the Rscript 01_annotate_common_variant_analysis_v2_07_06_2023.R for annotation.

 # plink \
 # --allow-no-sex \
 # --bfile merged_ccss \
 # --maf 0.01 \
 # --hwe 1e-06 \
 # --logistic hide-covar \
 # --extract assoc_test_two_vars.txt \
 # --ci 0.95 \
 # --pheno pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 # --pheno-name CMP \
 # --covar pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno \
 # --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 # --out sccss_org_ccss_exp_vars