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


## Annotate this plink file merged_plink_PLP
module load gatk/4.1.8.0
module load vt
module load vcftools
module load bcftools
module load tabix
module load zlib/1.2.5
module load java/13.0.1
module load vep/v108
module load samtools

export THREADS=2; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/annotation/snpeff/"; \
# export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/12_10_2023/clinvar.vcf.gz"; \
# export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
# export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/dbnsfp/new_dbNSFP4.4a_variant.chr${CHR}.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl";


plink --bfile merged_plink_PLP --keep-allele-order --recode vcf --out merged_plink_PLP_VCF
awk 'BEGIN {OFS="\t"} {if ($1 ~ /^#/){print $0}else{$1="chr"$1; print $0}}' merged_plink_PLP_VCF.vcf| less -S > merged_plink_PLP_VCF2.vcf

# awk 'BEGIN {OFS="\t"} {if ($1 ~ /^##contig/){gsub(/<ID=/,"<ID=chr"); print} else {print}}' merged_plink_PLP_VCF2.vcf > merged_plink_PLP_VCF.vcf

# In VCF file, replace 
##contig=<ID=1,length=201377675>
##contig=<ID=2,length=178807417>
##contig=<ID=3,length=52454023>
##contig=<ID=6,length=7586637>
##contig=<ID=10,length=119677818>
##contig=<ID=14,length=23435654>
##contig=<ID=17,length=39666549>

# with

##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>



module load htslib
bgzip -c merged_plink_PLP_VCF2.vcf > merged_plink_PLP_VCF2.vcf.gz
tabix -p vcf merged_plink_PLP_VCF2.vcf.gz

gatk VariantAnnotator \
   -R /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
   -V merged_plink_PLP_VCF2.vcf.gz \
   -L merged_plink_PLP_VCF2.vcf.gz \
   -D /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz \
   -O  merged_plink_PLP_VCF_ann_dbsnp.vcf

cat merged_plink_PLP_VCF_ann_dbsnp.vcf| grep -v "^##" | cut -f1-8 | awk -F'\t' '!_[$3]++' >> merged_plink_PLP_VCF_ann_dbsnp2.vcf
sed -i 's/\t\*\t/\t<*:DEL>\t/g' merged_plink_PLP_VCF_ann_dbsnp2.vcf
## Adding snpEff annotation; these are genome annotations from ENSEMBL, created from GRCh38/hg38 reference genome sequence
module load java/openjdk-11
${JAVA} ${JAVAOPTS} -jar ${SNPEFF} -v GRCh38.105  merged_plink_PLP_VCF_ann_dbsnp2.vcf > merged_plink_PLP_VCF_ann_dbsnp-snpeff.vcf


${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} merged_plink_PLP_VCF_ann_dbsnp-snpeff.vcf > merged_plink_PLP_VCF_ann_dbsnp-snpeff-ExAC.0.3.GRCh38.vcf

module load java/13.0.1
module load java/openjdk-11
# ## GnomAD with Exome
# ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/gnomAD/gnomad.genomes.v4.0.sites.${CHR}.vcf.bgz \
# merged_plink_PLP_VCF_ann_dbsnp-snpeff-ExAC.0.3.GRCh38.vcf > merged_plink_PLP_VCF_ann_dbsnp-snpeff-ExAC.0.3.GRCh38.gnomAD.vcf

# Input VCF file
INPUT_VCF="merged_plink_PLP_VCF_ann_dbsnp-snpeff-ExAC.0.3.GRCh38.vcf"
# Temporary output file
grep ^## ${INPUT_VCF} > HEADER
# Loop through chromosomes 1-22, X, Y
for CHR in 1 2 3 6 10 14 17; do
    # Define gnomAD file path for the current chromosome
    GNOMAD_VCF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/gnomAD/gnomad.genomes.v4.0.sites.chr${CHR}.vcf.bgz"
    grep ^chr${CHR} ${INPUT_VCF} > tmp_chr${CHR}_input.vcf
    cat HEADER tmp_chr${CHR}_input.vcf > chr${CHR}_input.vcf
    # Run SnpSift annotation and append to temporary output file
    ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate ${GNOMAD_VCF} chr${CHR}_input.vcf > chr${CHR}_out_gnomad.vcf
    DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/dbnsfp/new_dbNSFP4.4a_variant.chr${CHR}.gz"; 
    ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -n -v chr${CHR}_out_gnomad.vcf > chr${CHR}_out_gnomad-snpeff-dbnsfp.vcf
    ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} chr${CHR}_out_gnomad-snpeff-dbnsfp.vcf > chr${CHR}_out_gnomad-snpeff-dbnsfp.clinvar_12_10_2023.vcf
    cat chr${CHR}_out_gnomad-snpeff-dbnsfp.clinvar_12_10_2023.vcf |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "dbNSFP_gnomAD_genomes_AF" "dbNSFP_gnomAD_genomes_NFE_AF" "dbNSFP_gnomAD_genomes_AMI_AF" "dbNSFP_gnomAD_genomes_ASJ_AF" "dbNSFP_gnomAD_genomes_MID_AF" "dbNSFP_gnomAD_genomes_FIN_AF" "dbNSFP_gnomAD_genomes_SAS_AF" "dbNSFP_gnomAD_genomes_EAS_AF" "dbNSFP_gnomAD_genomes_AFR_AF" "dbNSFP_gnomAD_genomes_POPMAX_AF" "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_exomes_NFE_AF" "dbNSFP_gnomAD_exomes_ASJ_AF" "dbNSFP_gnomAD_exomes_FIN_AF" "dbNSFP_gnomAD_exomes_SAS_AF" "dbNSFP_gnomAD_exomes_EAS_AF" "dbNSFP_gnomAD_exomes_AFR_AF" "dbNSFP_gnomAD_exomes_POPMAX_AF" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_fin" "AF_raw" "AF" "AF_grpmax" "AF_joint_nfe" "AF_joint_afr" "AF_joint_eas" "AF_joint_sas" "AF_joint_fin" "AF_joint_raw" "AF_joint" "AF_grpmax_joint" > chr${CHR}_out_gnomad-snpeff-dbnsfp.clinvar_12_10_2023_FIELDS-simple.txt
done

## merge back
head -1 chr1_out_gnomad-snpeff-dbnsfp.clinvar_12_10_2023_FIELDS-simple.txt > HEADER1
for CHR in 1 2 3 6 10 14 17; do
tail -n +2 chr${CHR}_out_gnomad-snpeff-dbnsfp.clinvar_12_10_2023_FIELDS-simple.txt >> tmp_merged_ann.txt
done

cat HEADER1 tmp_merged_ann.txt > nine_genes_merged_snpeff_ann.txt

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


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_exp_samples.txt .


## Merge plink files
ls *all_vars.bim| grep bim| sed 's/.bim//' | sed -n '1d;p' > wanted_plink
plink --bfile ccss_exp_BAG3_all_vars --merge-list wanted_plink --keep-allele-order --out merged_plink_PLP

# subset to SJLIFE and CCSS_Exp EUR cohort
plink --bfile merged_plink_PLP --keep sjlife_ccss_exp_samples.txt --extract snpEFF_AF_eur_0.0001_vars.txt  --make-bed --keep-allele-order --out sjlife_ccss_exp_merged_EUR
plink --bfile merged_plink_PLP_sjlife --keep ../afr_samples --extract snpEFF_AF_afr_0.0001_vars.txt --make-bed --keep-allele-order --out sjlife_AFR

plink --bfile sjlife_ccss_exp_merged_EUR --freq -out sjlife_ccss_exp_merged_EUR.freq
plink --bfile sjlife_AFR --freq -out sjlife_AFR.freq
## Get allele count for each P/LP variants in nine gene in CCSS and SJLIFE using Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\rare_variant_analysis_v2\0.Annotation_v2_allele_count.R
## Now read these two files in Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\combined_cohorts\rare_variant_analysis\0.Annotation_v2; part 2


## Merge all SJLIFE maf 0.01 files and CCSS SJLIFE files and create raw files


## Extract ccss and sjlife overlapping variants 

## Run 0.Annotation_v2.R code to get the list of variants for analysis
plink --bfile sjlife_ccss_exp_merged_EUR --recodeA --out sjlife_ccss_exp_SNPS_maf_lt_0.0001_gnomad_recodeA_EUR
# snpEFF_EUR_vars_final_0.001_vars.txt is the same number of variants in EUR
plink --bfile sjlife_AFR --recodeA --out sjlife_SNPS_maf_lt_0.0001_gnomad_recodeA_AFR
# plink --bfile sjlife_AFR --geno 0.01 --extract snpEFF_EUR_vars_final_0.0001_vars.txt --recodeA --out sjlife_SNPS_maf_lt_0.0001_gnomad_recodeA_AFR


###############################
## Create PLP variant variables
# extract variants by gene categories
###############################
## Extract variants for P/LP 
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined

## EUR
table(AF_EUR.0.0001.final$`ANN[*].GENE`)
# BAG3   DSP  LMNA  MYH7 SCN5A  TCAP TNNT2   TTN 
# 10    21     9     5    18     4     2   320 

# length(TTN_PSI)
# # 267
# length(TTN_PSI_A_Band)
# # 157


cat BAG3_EUR_vars.txt DSP_EUR_vars.txt LMNA_EUR_vars.txt MYH7_EUR_vars.txt SCN5A_EUR_vars.txt TCAP_EUR_vars.txt TNNT2_EUR_vars.txt TTN_PSI_EUR_vars.txt > ALL_PLP_With_TTN_PSI_EUR.txt
cat BAG3_EUR_vars.txt DSP_EUR_vars.txt LMNA_EUR_vars.txt MYH7_EUR_vars.txt SCN5A_EUR_vars.txt TCAP_EUR_vars.txt TNNT2_EUR_vars.txt TTN_PSI_A_Band_EUR_vars.txt > ALL_PLP_With_TTN_PSI_A_Band_EUR.txt


plink --bfile sjlife_ccss_exp_merged_EUR --extract BAG3_EUR_vars.txt --out BAG3_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract DSP_EUR_vars.txt --out DSP_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract LMNA_EUR_vars.txt --out LMNA_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract MYH7_EUR_vars.txt --out MYH7_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract SCN5A_EUR_vars.txt --out SCN5A_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TCAP_EUR_vars.txt --out TCAP_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TNNT2_EUR_vars.txt --out TNNT2_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TTN_PSI_EUR_vars.txt --out TTN_PSI_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TTN_PSI_A_Band_EUR_vars.txt --out TTN_PSI_A_Band_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract BAG3.TTN_PSI_EUR_vars.txt --out BAG3.TTN_PSI_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract BAG3.TTN_PSI_A_Band_EUR_vars.txt --out BAG3.TTN_PSI_A_Band_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract ALL_PLP_With_TTN_PSI_EUR.txt --out ALL_PLP_With_TTN_PSI_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract ALL_PLP_With_TTN_PSI_A_Band_EUR.txt --out ALL_PLP_With_TTN_PSI_A_Band_EUR_vars_recodeA --recode A 





## AFR
# BAG3   DSP  LMNA  MYH7 SCN5A  TCAP TNNT2   TTN 
# 9    22     8     7    14     4     6   351 

# length(TTN_PSI)
# # 292
# length(TTN_PSI_A_Band)
# # 165


cat BAG3_AFR_vars.txt DSP_AFR_vars.txt LMNA_AFR_vars.txt MYH7_AFR_vars.txt SCN5A_AFR_vars.txt TCAP_AFR_vars.txt TNNT2_AFR_vars.txt TTN_PSI_AFR_vars.txt > ALL_PLP_With_TTN_PSI_AFR.txt
cat BAG3_AFR_vars.txt DSP_AFR_vars.txt LMNA_AFR_vars.txt MYH7_AFR_vars.txt SCN5A_AFR_vars.txt TCAP_AFR_vars.txt TNNT2_AFR_vars.txt TTN_PSI_A_Band_AFR_vars.txt > ALL_PLP_With_TTN_PSI_A_Band_AFR.txt


plink --bfile sjlife_AFR --extract BAG3_AFR_vars.txt --out BAG3_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract DSP_AFR_vars.txt --out DSP_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract LMNA_AFR_vars.txt --out LMNA_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract MYH7_AFR_vars.txt --out MYH7_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract SCN5A_AFR_vars.txt --out SCN5A_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TCAP_AFR_vars.txt --out TCAP_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TNNT2_AFR_vars.txt --out TNNT2_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TTN_PSI_AFR_vars.txt --out TTN_PSI_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TTN_PSI_A_Band_AFR_vars.txt --out TTN_PSI_A_Band_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract BAG3.TTN_PSI_AFR_vars.txt --out BAG3.TTN_PSI_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract BAG3.TTN_PSI_A_Band_AFR_vars.txt --out BAG3.TTN_PSI_A_Band_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract ALL_PLP_With_TTN_PSI_AFR.txt --out ALL_PLP_With_TTN_PSI_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract ALL_PLP_With_TTN_PSI_A_Band_AFR.txt --out ALL_PLP_With_TTN_PSI_A_Band_AFR_vars_recodeA --recode A 




#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################


# 10/10/2024
#############################
## Extract from PreQC data ##
#############################

# 1. Protein-Altering Variants (BAG3, LMNA, TCAP, TNNC1, TNNT2)
# These include missense, nonsense, and frameshift variants:
# Missense variants: missense_variant
# Nonsense variants: stop_gained
# Frameshift variants: frameshift_variant
# In-frame insertions/deletions: inframe_insertion, inframe_deletion

# 2. Missense Variants and In-frame Insertion/Deletion (MYH7)
# Missense variants: missense_variant
# In-frame variants: inframe_insertion, inframe_deletion

# 3. Frameshift, Stop-Gained, Splice-Donor, and Splice-Acceptor Variants (DSP, SCN5A, TTN)
# Frameshift variants: frameshift_variant
# Stop-gained variants: stop_gained
# Splice-donor/acceptor variants: splice_donor_variant, splice_acceptor_variant



# mkdir /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes_SNPEFF
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes_SNPEFF
# DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/"
# head -1 ${DIR}/chr2.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > header

# ## Extract variants regrions from nine genes: BAG3, DSP, LMNA, MYH7, SCN5A, TCAP, TNNC1, TNNT2, and TTN 
# GENE=NINE_GENES_ANNOVAR
# cat header > ${GENE}


# # cat ../ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| awk '/^chr10/ && ($2 >= 119651380 && $2 <= 119677819) && ($38+0 < 0.01 && $38 !=".") && ($44+0 < 0.01 && $44 !=".")' >> ${GENE}
# # BAG3; 10q26.11
# cat ${DIR}/chr10.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt| 
# awk '/^chr10/ && $2 >= 119651380 && $2 <= 119677819' >> ${GENE}
# # DSP
# cat ${DIR}/chr6.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
# awk '/^chr6/ && $2 >= 7541671 && $2 <= 7586714' >> ${GENE}
# # LMNA
# cat ${DIR}/chr1.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
# awk '/^chr1/ && $2 >= 156082573 && $2 <= 156140081' >> ${GENE}
# # MYH7
# cat ${DIR}/chr14.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
# awk '/^chr14/ && $2 >= 23412740 && $2 <= 23435660' >> ${GENE}
# # SCN5A
# cat ${DIR}/chr3.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
# awk '/^chr3/ && $2 >= 38548062 && $2 <= 38649687' >> ${GENE}
# # TCAP
# cat ${DIR}/chr17.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
# awk '/^chr17/ && $2 >= 39665349 && $2 <= 39666554' >> ${GENE}
# # TNNC1
# cat ${DIR}/chr3.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
# awk '/^chr3/ && $2 >= 52451100 && $2 <= 52454041' >> ${GENE}
# # TNNT2
# cat ${DIR}/chr1.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
# awk '/^chr1/ && $2 >= 201359014 && $2 <= 201377680' >> ${GENE}
# # TTN; 2q31.2
# cat ${DIR}/chr2.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt|
# awk '/^chr2/ && $2 >= 178525989 && $2 <= 178807423' >> ${GENE}


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined
# ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_exp_samples.txt .

# ## Merge plink files
# ls *all_vars.bim| grep bim| sed 's/.bim//' | sed -n '1d;p' > wanted_plink
# plink --bfile ccss_exp_BAG3_all_vars --merge-list wanted_plink --keep-allele-order --out merged_plink_PLP

# subset to SJLIFE and CCSS_Exp EUR cohort
plink --bfile merged_plink_PLP --keep sjlife_ccss_exp_samples.txt --extract snpEFF_AF_eur_0.0001_vars.txt  --make-bed --keep-allele-order --out sjlife_ccss_exp_merged_EUR
plink --bfile merged_plink_PLP_sjlife --keep ../afr_samples --extract snpEFF_AF_afr_0.0001_vars.txt --make-bed --keep-allele-order --out sjlife_AFR

plink --bfile sjlife_ccss_exp_merged_EUR --freq -out sjlife_ccss_exp_merged_EUR.freq
plink --bfile sjlife_AFR --freq -out sjlife_AFR.freq
## Get allele count for each P/LP variants in nine gene in CCSS and SJLIFE using Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\rare_variant_analysis_v2\0.Annotation_v2_allele_count.R
## Now read these two files in Rscipt Z:\ResearchHome\ClusterHome\aneupane\St_Jude\Achal_St_Jude\rcodes\TITN_BAG3\combined_cohorts\rare_variant_analysis\0.Annotation_v2; part 2


## Merge all SJLIFE maf 0.01 files and CCSS SJLIFE files and create raw files


## Extract ccss and sjlife overlapping variants 

## Run 0.Annotation_v2.R code to get the list of variants for analysis
plink --bfile sjlife_ccss_exp_merged_EUR --recodeA --out sjlife_ccss_exp_SNPS_maf_lt_0.0001_gnomad_recodeA_EUR
# snpEFF_EUR_vars_final_0.001_vars.txt is the same number of variants in EUR
plink --bfile sjlife_AFR --recodeA --out sjlife_SNPS_maf_lt_0.0001_gnomad_recodeA_AFR
# plink --bfile sjlife_AFR --geno 0.01 --extract snpEFF_EUR_vars_final_0.0001_vars.txt --recodeA --out sjlife_SNPS_maf_lt_0.0001_gnomad_recodeA_AFR


###############################
## Create PLP variant variables
# extract variants by gene categories
###############################
## Extract variants for P/LP 
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined

## EUR
table(AF_EUR.0.0001.final$`ANN[*].GENE`)



cat BAG3_EUR_vars.txt DSP_EUR_vars.txt LMNA_EUR_vars.txt MYH7_EUR_vars.txt SCN5A_EUR_vars.txt TCAP_EUR_vars.txt TNNC1_EUR_vars.txt TNNT2_EUR_vars.txt TTN_PSI_EUR_vars.txt > ALL_PLP_With_TTN_PSI_EUR.txt
cat BAG3_EUR_vars.txt DSP_EUR_vars.txt LMNA_EUR_vars.txt MYH7_EUR_vars.txt SCN5A_EUR_vars.txt TCAP_EUR_vars.txt TNNC1_EUR_vars.txt TNNT2_EUR_vars.txt TTN_PSI_A_Band_EUR_vars.txt > ALL_PLP_With_TTN_PSI_A_Band_EUR.txt


plink --bfile sjlife_ccss_exp_merged_EUR --extract BAG3_EUR_vars.txt --out BAG3_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract DSP_EUR_vars.txt --out DSP_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract LMNA_EUR_vars.txt --out LMNA_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract MYH7_EUR_vars.txt --out MYH7_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract SCN5A_EUR_vars.txt --out SCN5A_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TCAP_EUR_vars.txt --out TCAP_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TNNC1_EUR_vars.txt --out TNNC1_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TNNT2_EUR_vars.txt --out TNNT2_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TTN_PSI_EUR_vars.txt --out TTN_PSI_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract TTN_PSI_A_Band_EUR_vars.txt --out TTN_PSI_A_Band_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract BAG3.TTN_PSI_EUR_vars.txt --out BAG3.TTN_PSI_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract BAG3.TTN_PSI_A_Band_EUR_vars.txt --out BAG3.TTN_PSI_A_Band_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract ALL_PLP_With_TTN_PSI_EUR.txt --out ALL_PLP_With_TTN_PSI_EUR_vars_recodeA --recode A 
plink --bfile sjlife_ccss_exp_merged_EUR --extract ALL_PLP_With_TTN_PSI_A_Band_EUR.txt --out ALL_PLP_With_TTN_PSI_A_Band_EUR_vars_recodeA --recode A 




cat BAG3_AFR_vars.txt DSP_AFR_vars.txt LMNA_AFR_vars.txt MYH7_AFR_vars.txt SCN5A_AFR_vars.txt TCAP_AFR_vars.txt TNNC1_AFR_vars.txt TNNT2_AFR_vars.txt TTN_PSI_AFR_vars.txt > ALL_PLP_With_TTN_PSI_AFR.txt
cat BAG3_AFR_vars.txt DSP_AFR_vars.txt LMNA_AFR_vars.txt MYH7_AFR_vars.txt SCN5A_AFR_vars.txt TCAP_AFR_vars.txt TNNC1_AFR_vars.txt TNNT2_AFR_vars.txt TTN_PSI_A_Band_AFR_vars.txt > ALL_PLP_With_TTN_PSI_A_Band_AFR.txt


plink --bfile sjlife_AFR --extract BAG3_AFR_vars.txt --out BAG3_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract DSP_AFR_vars.txt --out DSP_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract LMNA_AFR_vars.txt --out LMNA_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract MYH7_AFR_vars.txt --out MYH7_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract SCN5A_AFR_vars.txt --out SCN5A_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TCAP_AFR_vars.txt --out TCAP_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TNNC1_AFR_vars.txt --out TNNC1_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TNNT2_AFR_vars.txt --out TNNT2_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TTN_PSI_AFR_vars.txt --out TTN_PSI_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract TTN_PSI_A_Band_AFR_vars.txt --out TTN_PSI_A_Band_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract BAG3.TTN_PSI_AFR_vars.txt --out BAG3.TTN_PSI_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract BAG3.TTN_PSI_A_Band_AFR_vars.txt --out BAG3.TTN_PSI_A_Band_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract ALL_PLP_With_TTN_PSI_AFR.txt --out ALL_PLP_With_TTN_PSI_AFR_vars_recodeA --recode A 
plink --bfile sjlife_AFR --extract ALL_PLP_With_TTN_PSI_A_Band_AFR.txt --out ALL_PLP_With_TTN_PSI_A_Band_AFR_vars_recodeA --recode A 



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

