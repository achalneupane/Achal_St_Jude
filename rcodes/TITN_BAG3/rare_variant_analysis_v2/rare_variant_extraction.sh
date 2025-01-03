cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes
DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar"
head -1 ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > header

## Extract variants regrions from nine genes: BAG3, DSP, LMNA, MYH7, SCN5A, TCAP, TNNC1, TNNT2, and TTN 
GENE=NINE_GENES_ANNOVAR
cat header > ${GENE}


# cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| awk '/^chr10/ && ($2 >= 119651380 && $2 <= 119677819) && ($38+0 < 0.01 && $38 !=".") && ($44+0 < 0.01 && $44 !=".")' >> ${GENE}
# BAG3; 10q26.11
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| 
awk '/^chr10/ && $2 >= 119651380 && $2 <= 119677819' >> ${GENE}
# DSP
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr6.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr6/ && $2 >= 7541671 && $2 <= 7586714' >> ${GENE}
# LMNA
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr1/ && $2 >= 156082573 && $2 <= 156140081' >> ${GENE}
# MYH7
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr14.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr14/ && $2 >= 23412740 && $2 <= 23435660' >> ${GENE}
# SCN5A
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr3.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr3/ && $2 >= 38548062 && $2 <= 38649687' >> ${GENE}
# TCAP
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr17.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr17/ && $2 >= 39665349 && $2 <= 39666554' >> ${GENE}
# TNNC1
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr3.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr3/ && $2 >= 52451100 && $2 <= 52454041' >> ${GENE}
# TNNT2
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr1/ && $2 >= 201359014 && $2 <= 201377680' >> ${GENE}
# TTN; 2q31.2
cat ${DIR}/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr2/ && $2 >= 178525989 && $2 <= 178807423' >> ${GENE}


##############################################
## Extract these nine genes from plink data ##
##############################################
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
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
--max-maf 0.01 \
--out sjlife_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz

plink \
--chr ${CHR} \
--from-bp ${START} \
--make-bed \
--max-maf 0.01 \
--out ccss_exp_${GENE}_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz



cut -d' ' -f1-2 /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3.pheno > sjlife_samples
cut -d' ' -f1-2 /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_exp_eur_cardiotoxic_exposed.pheno > ccss_exp_samples


plink1.9 --bfile ${BFILE} --keep amyloid_imaging_sample_list.txt --make-bed --keep-allele-order --out ${BFILE}_Amyloid_Imaging


# subset to SJLIFE EUR cohort
for file in ls sjlife*0.01.bim; do
BFILE="$(echo ${file%.*})"
plink --bfile ${BFILE} --keep sjlife_samples --make-bed --keep-allele-order --max-maf 0.01 --out ${BFILE}_final
done

# subset to CCSS EUR cohort
for file in ls ccss_exp*0.01.bim; do
BFILE="$(echo ${file%.*})"
plink --bfile ${BFILE} --keep ccss_exp_samples --make-bed --keep-allele-order --max-maf 0.01 --out ${BFILE}_final
done

cat ccss_exp*0.01_final.bim > all_ccss_exp_vars_lt_maf_0.01.txt
cat sjlife*0.01_final.bim > all_sjlife_vars_lt_maf_0.01.txt

## Get allele count for each P/LP variants in nine gene in CCSS and SJLIFE using Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\rare_variant_analysis_v2\0.Annotation_v2_allele_count.R
## Now read these two files in Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\rare_variant_analysis_v2\0.Annotation_v2; part 2


## Merge all SJLIFE maf 0.01 files and CCSS SJLIFE files and create raw files

##################################################
## Extract ccss and sjlife overlapping variants ##
##################################################

ls sjlife*0.01_final.bim| grep bim| sed 's/.bim//' | sed -n '1d;p' > sjlife_wanted_plink
plink --bfile sjlife_BAG3_maf_lt_0.01_final --merge-list sjlife_wanted_plink --out sjlife_all_genes_maf_lt_0.01_final
plink --bfile sjlife_all_genes_maf_lt_0.01_final --extract sjlife_SNPS_maf_lt_0.01_gnomad_also_common_in_ccss.txt --recodeA --out sjlife_SNPS_maf_lt_0.01_gnomad_also_common_in_ccss_recodeA

ls ccss*0.01_final.bim| grep bim| sed 's/.bim//' | sed -n '1d;p' > ccss_wanted_plink

## Before
plink --bfile ccss_exp_BAG3_maf_lt_0.01_final --merge-list ccss_wanted_plink --out ccss_all_genes_maf_lt_0.01_final
plink --bfile ccss_all_genes_maf_lt_0.01_final --extract ccss_SNPS_maf_lt_0.01_gnomad_also_common_in_sjlife.txt --recodeA --out ccss_SNPS_maf_lt_0.01_gnomad_also_common_in_sjlife_recodeA


## Now, after removing the younger samples within 5 years of primary cancer
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes

plink --bfile ccss_all_genes_maf_lt_0.01_final --remove /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_exp.young.CA.to.remove.txt --make-bed --out ccss_all_genes_maf_lt_0.01_final_without_younger_samples
plink --bfile ccss_all_genes_maf_lt_0.01_final_without_younger_samples --extract ccss_SNPS_maf_lt_0.01_gnomad_also_common_in_sjlife.txt --recodeA --out ccss_SNPS_maf_lt_0.01_gnomad_also_common_in_sjlife_without_younger_samples_recodeA






## Rare P/LP analysis for SJLIFE (EUR and African) echo 
# First merge all rare PLP plink files for SJLIFE together.
echo "sjlife_BAG3_maf_lt_0.01" > merge_list.txt
echo "sjlife_DSP_maf_lt_0.01" >> merge_list.txt
echo "sjlife_LMNA_maf_lt_0.01" >> merge_list.txt
echo "sjlife_MYH7_maf_lt_0.01" >> merge_list.txt
echo "sjlife_SCN5A_maf_lt_0.01" >> merge_list.txt
echo "sjlife_TCAP_maf_lt_0.01" >> merge_list.txt
echo "sjlife_TNNC1_maf_lt_0.01" >> merge_list.txt
echo "sjlife_TNNT2_maf_lt_0.01" >> merge_list.txt
echo "sjlife_TTN_maf_lt_0.01" >> merge_list.txt

plink --merge-list merge_list.txt --make-bed --extract sjlife_SNPS_maf_lt_0.01_gnomad_also_common_in_ccss.txt --out merged_data_maf_lt_0.01
## Keep with carriers only
plink --bfile merged_data_maf_lt_0.01 --out ALL_sjlife_samples_SNPS_maf_lt_0.01_gnomad_also_common_in_ccss_recodeA --recode A
# LMNA
chr1:156138821:C:T 
# TNTT2
chr1:201361214:C:T
chr1:201364335:CG:C
# TTN
chr2:178542266:A:G
chr2:178577602:C:T
chr2:178615321:A:G
chr2:178630241:G:A
chr2:178651534:G:T
chr2:178664926:G:C
chr2:178665777:G:A
chr2:178669675:G:T
chr2:178692016:C:T
chr2:178693609:C:T
chr2:178701528:C:G
chr2:178702065:T:A
chr2:178728776:G:T
chr2:178746079:G:A
chr2:178782806:C:T
# BAG3
chr10:119672255:C:T

(base) [aneupane@splprhpc11 pablo_garcia_et_al_nine_genes]$ vi LMNA_sjlife_vars
(base) [aneupane@splprhpc11 pablo_garcia_et_al_nine_genes]$ vi TNTT2_sjlife_vars
(base) [aneupane@splprhpc11 pablo_garcia_et_al_nine_genes]$ vi TTN_sjlife_vars
(base) [aneupane@splprhpc11 pablo_garcia_et_al_nine_genes]$ vi BAG3_sjlife_vars
(base) [aneupane@splprhpc11 pablo_garcia_et_al_nine_genes]$ vi ALL_PLP_sjlife_vars

plink --bfile merged_data_maf_lt_0.01 --extract LMNA_sjlife_vars --out LMNA_PLP_SJLIFE_recodeA --recode A
plink --bfile merged_data_maf_lt_0.01 --extract TNTT2_sjlife_vars --out TNTT2_PLP_SJLIFE_recodeA --recode A
plink --bfile merged_data_maf_lt_0.01 --extract TTN_sjlife_vars --out TTN_PLP_SJLIFE_recodeA --recode A
plink --bfile merged_data_maf_lt_0.01 --extract BAG3_sjlife_vars --out BAG3_PLP_SJLIFE_recodeA --recode A
plink --bfile merged_data_maf_lt_0.01 --extract ALL_PLP_sjlife_vars --out ALL_PLP_SJLIFE_recodeA --recode A



# SJL1751413
# SJL4172706
# SJL5203001

# SJL1727210 # CTCAE grade 0