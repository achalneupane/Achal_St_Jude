# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes
# DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar"
# head -1 ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > header

# ## Extract variants regrions from nine genes: BAG3, DSP, LMNA, MYH7, SCN5A, TCAP, TNNC1, TNNT2, and TTN 
# GENE=NINE_GENES_ANNOVAR
# cat header > ${GENE}


# # cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| awk '/^chr10/ && ($2 >= 119651380 && $2 <= 119677819) && ($38+0 < 0.01 && $38 !=".") && ($44+0 < 0.01 && $44 !=".")' >> ${GENE}
# # BAG3; 10q26.11
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| 
# awk '/^chr10/ && $2 >= 119651380 && $2 <= 119677819' >> ${GENE}
# # DSP
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr6.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
# awk '/^chr6/ && $2 >= 7541671 && $2 <= 7586714' >> ${GENE}
# # LMNA
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
# awk '/^chr1/ && $2 >= 156082573 && $2 <= 156140081' >> ${GENE}
# # MYH7
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr14.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
# awk '/^chr14/ && $2 >= 23412740 && $2 <= 23435660' >> ${GENE}
# # SCN5A
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr3.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
# awk '/^chr3/ && $2 >= 38548062 && $2 <= 38649687' >> ${GENE}
# # TCAP
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr17.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
# awk '/^chr17/ && $2 >= 39665349 && $2 <= 39666554' >> ${GENE}
# # TNNC1
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr3.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
# awk '/^chr3/ && $2 >= 52451100 && $2 <= 52454041' >> ${GENE}
# # TNNT2
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
# awk '/^chr1/ && $2 >= 201359014 && $2 <= 201377680' >> ${GENE}
# # TTN; 2q31.2
# cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
# awk '/^chr2/ && $2 >= 178525989 && $2 <= 178807423' >> ${GENE}


##############################################
## Extract these nine genes from plink data ##
##############################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined
module load plink/1.90b
## BAG3
CHR=10
GENE=BAG3
START=119651380
END=119677819
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz

## DSP
CHR=6
GENE=DSP
START=7541671
END=7586714
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz


## LMNA
CHR=1
GENE=LMNA
START=156082573
END=156140081
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz


## MYH7
CHR=14
GENE=MYH7
START=23412740
END=23435660
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz


## SCN5A
CHR=3
GENE=SCN5A
START=38548062
END=38649687
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz


## TCAP
CHR=17
GENE=TCAP
START=39665349
END=39666554
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz


## TNNC1
CHR=3
GENE=TNNC1
START=52451100
END=52454041
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz


## TNNT2
CHR=1
GENE=TNNT2
START=201359014
END=201377680
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz


## TTN
CHR=2
GENE=TTN
START=178525989
END=178807423
plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out sjlife_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--out ccss_exp_${GENE}_all_vars \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz



ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_exp_samples.txt .

## Merge plink files
ls *all_vars.bim| grep bim| sed 's/.bim//' | sed -n '1d;p'| grep sjlife > wanted_plink
plink --bfile sjlife_BAG3_all_vars --merge-list wanted_plink_sjlife --keep-allele-order --out merged_plink_PLP_sjlife

ls *all_vars.bim| grep bim| sed 's/.bim//' | sed -n '1d;p'| grep ccss > wanted_plink
plink --bfile ccss_exp_BAG3_all_vars --merge-list wanted_plink_ccss_exp --keep-allele-order --out merged_plink_PLP_ccss_exp

## Check for allele flips
ls *all_vars.bim| grep bim| sed 's/.bim//' | sed -n '1d;p' > wanted_plink
plink --bfile ccss_exp_BAG3_all_vars --merge-list wanted_plink --keep-allele-order --out merged_plink_PLP

## Get allele count for each P/LP variants in nine gene in CCSS and SJLIFE using Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\rare_variant_analysis_v2\0.Annotation_v2_allele_count.R
## Now read these two files in Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\combined_cohorts\rare_variant_analysis\0.Annotation_v2; part 2


## Merge all SJLIFE maf 0.01 files and CCSS SJLIFE files and create raw files

##################################################################################################################




# 09/16/2024
#############################
## Extract from PreQC data ##
#############################
mkdir /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes_SNPEFF
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes_SNPEFF
DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/"
head -1 ${DIR}/chr2.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > header

## Extract variants regrions from nine genes: BAG3, DSP, LMNA, MYH7, SCN5A, TCAP, TNNC1, TNNT2, and TTN 
GENE=NINE_GENES_ANNOVAR
cat header > ${GENE}


# cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| awk '/^chr10/ && ($2 >= 119651380 && $2 <= 119677819) && ($38+0 < 0.01 && $38 !=".") && ($44+0 < 0.01 && $44 !=".")' >> ${GENE}
# BAG3; 10q26.11
cat ${DIR}/chr10.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt| 
awk '/^chr10/ && $2 >= 119651380 && $2 <= 119677819' >> ${GENE}
# DSP
cat ${DIR}/chr6.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
awk '/^chr6/ && $2 >= 7541671 && $2 <= 7586714' >> ${GENE}
# LMNA
cat ${DIR}/chr1.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
awk '/^chr1/ && $2 >= 156082573 && $2 <= 156140081' >> ${GENE}
# MYH7
cat ${DIR}/chr14.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
awk '/^chr14/ && $2 >= 23412740 && $2 <= 23435660' >> ${GENE}
# SCN5A
cat ${DIR}/chr3.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
awk '/^chr3/ && $2 >= 38548062 && $2 <= 38649687' >> ${GENE}
# TCAP
cat ${DIR}/chr17.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
awk '/^chr17/ && $2 >= 39665349 && $2 <= 39666554' >> ${GENE}
# TNNC1
cat ${DIR}/chr3.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
awk '/^chr3/ && $2 >= 52451100 && $2 <= 52454041' >> ${GENE}
# TNNT2
cat ${DIR}/chr1.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
awk '/^chr1/ && $2 >= 201359014 && $2 <= 201377680' >> ${GENE}
# TTN; 2q31.2
cat ${DIR}/chr2.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
awk '/^chr2/ && $2 >= 178525989 && $2 <= 178807423' >> ${GENE}


ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_exp_samples.txt .

## Merge plink files
ls *all_vars.bim| grep bim| sed 's/.bim//' | sed -n '1d;p' > wanted_plink
plink --bfile ccss_exp_BAG3_all_vars --merge-list wanted_plink --keep-allele-order --out merged_plink_PLP

# subset to SJLIFE and CCSS_Exp EUR cohort
plink --bfile merged_plink_PLP --keep sjlife_ccss_exp_samples.txt --extract snpEFF_AF_eur_0.001_vars.txt  --make-bed --keep-allele-order --out sjlife_ccss_exp_merged_EUR
plink --bfile merged_plink_PLP_sjlife --keep ../afr_samples --extract snpEFF_AF_eur_0.001_vars.txt --make-bed --keep-allele-order --out sjlife_AFR

## Get allele count for each P/LP variants in nine gene in CCSS and SJLIFE using Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\rare_variant_analysis_v2\0.Annotation_v2_allele_count.R
## Now read these two files in Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\combined_cohorts\rare_variant_analysis\0.Annotation_v2; part 2


## Merge all SJLIFE maf 0.01 files and CCSS SJLIFE files and create raw files


## Extract ccss and sjlife overlapping variants 

## Run 0.Annotation_v2.R code to get the list of variants for analysis
plink --bfile sjlife_ccss_exp_merged_EUR --geno 0.01 --recodeA --out sjlife_ccss_exp_SNPS_maf_lt_0.001_gnomad_recodeA_EUR
# snpEFF_EUR_vars_final_0.001_vars.txt is the same number of variants in EUR
plink --bfile sjlife_AFR --geno 0.01 --extract snpEFF_EUR_vars_final_0.001_vars.txt --recodeA --out sjlife_SNPS_maf_lt_0.001_gnomad_recodeA_AFR

# ## AFR specific variant
# # plink --bfile merged_plink_PLP_sjlife --keep ../afr_samples --extract snpEFF_AF_afr_0.001_vars.txt --make-bed --keep-allele-order --out sjlife_AFR
# plink --bfile sjlife_AFR --geno 0.01 --recodeA --out sjlife_SNPS_maf_lt_0.001_gnomad_recodeA_AFR





# ## Extract these nine genes from plink data from preQC data 
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined
# module load plink/1.90b
# ## BAG3
# CHR=10
# GENE=BAG3
# START=119651380
# END=119677819
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz

# ## DSP
# CHR=6
# GENE=DSP
# START=7541671
# END=7586714
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz


# ## LMNA
# CHR=1
# GENE=LMNA
# START=156082573
# END=156140081
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz


# ## MYH7
# CHR=14
# GENE=MYH7
# START=23412740
# END=23435660
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz


# ## SCN5A
# CHR=3
# GENE=SCN5A
# START=38548062
# END=38649687
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz


# ## TCAP
# CHR=17
# GENE=TCAP
# START=39665349
# END=39666554
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz


# ## TNNC1
# CHR=3
# GENE=TNNC1
# START=52451100
# END=52454041
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz


# ## TNNT2
# CHR=1
# GENE=TNNT2
# START=201359014
# END=201377680
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz


# ## TTN
# CHR=2
# GENE=TTN
# START=178525989
# END=178807423
# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out sjlife_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC.vcf.gz

# plink \
# --chr ${CHR} \
# --from-bp ${START} \
# --make-bed \
# --out ccss_exp_${GENE}_all_vars_preQC \
# --to-bp ${END} \
# --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr${CHR}_ID_edited.vcf.gz

