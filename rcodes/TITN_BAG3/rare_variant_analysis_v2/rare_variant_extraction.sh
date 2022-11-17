cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/pablo_garcia_et_al_nine_genes
head -1 ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > header

## Extract variants regrions from nine genes: BAG3, DSP, LMNA, MYH7, SCN5A, TCAP, TNNC1, TNNT2, and TTN 
GENE=NINE_GENES
cat header > ${GENE}
# cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| awk '/^chr10/ && ($2 >= 119651380 && $2 <= 119677819) && ($38+0 < 0.01 && $38 !=".") && ($44+0 < 0.01 && $44 !=".")' >> ${GENE}
# BAG3; 10q26.11
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| 
awk '/^chr10/ && $2 >= 119651380 && $2 <= 119677819' >> ${GENE}
# DSP
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr6.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr6/ && $2 >= 7541671 && $2 <= 7586714' >> ${GENE}
# LMNA
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr1/ && $2 >= 156082573 && $2 <= 156140081' >> ${GENE}
# MYH7
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr14.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr14/ && $2 >= 23412740 && $2 <= 23435660' >> ${GENE}
# SCN5A
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr3.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr3/ && $2 >= 38548062 && $2 <= 38649687' >> ${GENE}
# TCAP
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr17.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr17/ && $2 >= 39665349 && $2 <= 39666554' >> ${GENE}
# TNNC1
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr3.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr3/ && $2 >= 52451100 && $2 <= 52454041' >> ${GENE}
# TNNT2
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
awk '/^chr1/ && $2 >= 201359014 && $2 <= 201377680' >> ${GENE}
# TTN; 2q31.2
cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt|
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
--out ccss_exp_bag3_maf_lt_0.01 \
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
--out ccss_exp_bag3_maf_lt_0.01 \
--to-bp ${END} \
--vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz
