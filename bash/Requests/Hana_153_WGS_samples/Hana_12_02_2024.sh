## Email on 12/2/2024
Hi Achal,

Happy Monday! Hope you had a great holiday!

I had a question about the WGS for some of the samples we sent you before. Would it be possible to send me the list of variants for the following samples (not sure if 25-3 and 26-3 were included in the ones we sent you)?

19-3
25-3
26-3
GW30
GW64
GW159
GW10
GW53
GW168
GW129
GW167
GW169

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples
module load bcftools

## Fix multiallelic and normalize
bcftools norm -Ou -m -any CAB_5428.haplotype.vcf.gz \
 | bcftools norm -Ou -f /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
 | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' \
 | bcftools view -Oz \
 > CAB_5428.haplotype.recalibrated.decomposed.vcf.gz


for chr in {1..22} X Y;do
  export chr=${chr}
  bsub -q priority -P HanaWGS -J QC1.chr$chr -eo logs/QC1.$chr.err -oo logs/QC1.$chr.out -R "rusage[mem=20000]" bash process_WGS_data_perchr.sh $chr
done



#######################
#!/bin/bash

# Processes performed to QC the SJLIFE WGS data provided by Comp. Bio. ["SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714.vcf.gz"]

# Load software packages
module load vcftools/0.1.13
module load bcftools/1.17
#module load vt/2016.11.07
module load plink/1.90b
# module load R/3.4.0
module load R/4.2.2-rhel8
module load bedtools/2.25.0
module load htslib/1.3.1

# Define/create files/directories
INDIR=//research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/

OUTDIR=//research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/

VCFROOT="CAB_5428.haplotype.recalibrated.decomposed"

stats=//research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/stats/

cd $OUTDIR

# Process data per chromosome
chr=$1

echo "Doing chromosome chr${chr}"
# Split data for each chromosome
vcftools --gzvcf ${INDIR}/${VCFROOT}.vcf.gz --chr chr$chr --recode  --stdout | bgzip  > ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz
# bcftools view -S ^${OUTDIR}/chrALL.Survivor_WGS.GATK4180.hg38_missingness.imiss.0.05plus_renamed -r chr$chr -Oz -o ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz ${INDIR}/chrALL.${VCFROOT}.vcf.gz

# Index the VCF file
tabix -pvcf ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz

# Write stats
bcftools stats ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz > ${stats}/${VCFROOT}_chr${chr}.bcftools_stats

# Run basic QC - sequence level
vcftools --gzvcf ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz \
--chr chr$chr \
--keep-filtered PASS \
--minGQ 20 \
--min-meanDP 10 \
--recode \
--stdout \
| bgzip > ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz

# bcftools view -f PASS -i 'MIN(GQ)>20 & MEAN(DP)>10' -Oz -o ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz

# Write stats
bcftools stats ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz > ${stats}/${VCFROOT}_chr${chr}.PASS.bcftools_stats

tabix -pvcf ${INDIR}/${VCFROOT}_chr${chr}.PASS.gz

# Also write PLINK files for further processing
plink --vcf ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz \
 --keep-allele-order \
 --allow-extra-chr 0 \
 --hwe 1e-10 \
 --make-bed --out ${VCFROOT}_chr${chr}.PASS.qced

# Perform LD pruning to obtain independent set of markers
plink --nonfounders \
 --bfile ${VCFROOT}_chr${chr}.PASS.qced \
 --exclude range /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/high-LD-regions-GRCh38.txt \
 --geno 0.01 \
 --hwe 0.0001 \
 --maf 0.1 \
 --indep-pairwise 100 25 0.2 \
 --out ${VCFROOT}_chr${chr}.PASS.qced_common_pruned
plink --nonfounders \
 --bfile ${VCFROOT}_chr${chr}.PASS.qced \
 --extract ${VCFROOT}_chr${chr}.PASS.qced_common_pruned.prune.in \
 --make-bed \
 --out ${VCFROOT}_chr${chr}.PASS_common_pruned_indep


#######################


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/Hana_request_12_2_2024
module load plink/1.90b

sed -i 's/\r$//' gene_regions.txt 
# (base) [aneupane@splprhpc12 Hana_request_12_2_2024]$ awk 'NR>1 {print "chr"$2, $3, $4, $1}' gene_regions.txt| head -5
# chr17 7217125 7225266 ACADVL
# chr17 63477061 63498380 ACE
# chr15 34790107 34795589 ACTC1
# chr1 236664141 236764631 ACTN2
# chr10 114043866 114046904 ADRB1

awk 'NR>1 {print "chr"$2, $3, $4, $1}' gene_regions.txt > gene_regions.bed
cat -A gene_regions.bed | head

plink --vcf ../CAB_5428.haplotype.recalibrated.decomposed.vcf.gz \
      --extract range gene_regions.bed \
      --make-bed --out WGS_153_samples \
      --vcf-half-call m --keep-allele-order --double-id

sed -i 's/^23/X/' WGS_153_samples.bim
plink --bfile WGS_153_samples --recodeA --keep-allele-order --out WGS_153_samples_recodeA 

# plink --vcf ../CAB_5428.haplotype.recalibrated.decomposed_chr22.vcf.gz       --extract range gene_regions.bed       --make-bed --out WGS_153_samples       --vcf-half-call m --keep-allele-order --double-id
plink --bfile WGS_153_samples --recode vcf --out WGS_153_samples_converted
sed -i 's/^23/X/' WGS_153_samples_converted.vcf 

## Since gnomad data is per chromosome we need to split the VCF per chromosome
grep "^#" WGS_153_samples_converted-annot-snpeff-ExAC.0.3.GRCh38.vcf > header.vcf
for chr in {1..22} X; do
  grep -P "^#|\bchr${chr}\b" WGS_153_samples_converted.vcf  > chr${chr}.vcf
done

####################
## Run annotation ##
####################
for i in {1..22} X; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export INPUT_VCF="${CHR}.vcf"; \
export THREADS=2; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/Hana_request_12_2_2024/"; \
# export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export REF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/12_10_2023/clinvar.vcf.gz"; \
# export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
# export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/dbnsfp/new_dbNSFP4.4a_variant.${CHR}.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
 -P "${CHR}_annotate" \
 -q "standard" \
 -J "${CHR}_HANA" \
 -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
 -n "${THREADS}" \
 -R "rusage[mem=40GB]" \
 "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/Hana_request_12_2_2024/annotate.sh"; \
done;

#!/usr/bin/bash
module load gatk/4.1.8.0
module load vt
module load vcftools
module load bcftools
module load tabix
module load zlib/1.2.5
module load java/13.0.1
module load java/openjdk-11
module load vep/v108
module load samtools


cd ${WORKDIR}

VCF="${INPUT_VCF}"

MAX_HEADER_LINES=5000
ANNOT_SOURCE="${VCF}"; ANNOT_PROJECT="${ANNOT_SOURCE%.*}-annot"

## Start annotating
# zcat ${ANNOT_SOURCE} | head -${MAX_HEADER_LINES} | grep "^##" > ${ANNOT_PROJECT}.vcf
cat ${ANNOT_SOURCE}| grep -v "^##" | cut -f1-8 | awk -F'\t' '!_[$3]++' >> ${ANNOT_PROJECT}.vcf
sed -i 's/\t\*\t/\t<*:DEL>\t/g' ${ANNOT_PROJECT}.vcf

## Adding snpEff annotation; these are genome annotations from ENSEMBL, created from GRCh38/hg38 reference genome sequence
${JAVA} ${JAVAOPTS} -jar ${SNPEFF} -v GRCh38.105  ${ANNOT_PROJECT}.vcf > ${ANNOT_PROJECT}-snpeff.vcf
# Adding EXaC db
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf
ANNOTATED="${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf"


## GnomAD with Genome
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/gnomAD/genome_4.1/gnomad.genomes.v4.1.sites.${CHR}.vcf.bgz \
"${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf" > "${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.with.gnomAD.vcf"
ANNOTATED="${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.with.gnomAD.vcf"

## Revel
vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
--offline --everything --force_overwrite --vcf \
--fasta /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
-i ${CHR}-annot-snpeff-ExAC.0.3.GRCh38.with.gnomAD.vcf \
-o ${CHR}.WGS_VCF_revel.vcf \
--assembly GRCh38 \
--dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/VEP_plugins-release-110 \
--plugin REVEL,/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/new_tabbed_revel_grch38.tsv.gz


## Adding dbNSFP
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -n -v ${CHR}.WGS_VCF_revel.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
ANNOTATED=${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
## Adding clinvar CLNSIG
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOTATED} > "$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf"

ANNOTATED="$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf"
# grep ^## new_chr18.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp.vcf| grep -i dbnsfp| egrep -vi 'gnomad|_AF|_AC ' ''|grep -o 'dbNSFP_[^,]*'| sort | sed 's/.*/"&"/'| tr '\n'
# cat ${ANNOTATED}|${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF" "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_aaalt" "dbNSFP_aapos" "dbNSFP_aaref" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_raw" "AF" > ${WORKDIR}/"$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf-FIELDS-simple2.txt"
OUTFILE="${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple.txt"
cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "dbNSFP_gnomAD_genomes_AF" "dbNSFP_gnomAD_genomes_NFE_AF" "dbNSFP_gnomAD_genomes_AMI_AF" "dbNSFP_gnomAD_genomes_ASJ_AF" "dbNSFP_gnomAD_genomes_MID_AF" "dbNSFP_gnomAD_genomes_FIN_AF" "dbNSFP_gnomAD_genomes_SAS_AF" "dbNSFP_gnomAD_genomes_EAS_AF" "dbNSFP_gnomAD_genomes_AFR_AF" "dbNSFP_gnomAD_genomes_POPMAX_AF" "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_exomes_NFE_AF" "dbNSFP_gnomAD_exomes_ASJ_AF" "dbNSFP_gnomAD_exomes_FIN_AF" "dbNSFP_gnomAD_exomes_SAS_AF" "dbNSFP_gnomAD_exomes_EAS_AF" "dbNSFP_gnomAD_exomes_AFR_AF" "dbNSFP_gnomAD_exomes_POPMAX_AF" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_fin" "AF_raw" "AF" > ${WORKDIR}/${OUTFILE}

# combine all annotation files
awk 'NR == 1 {print; next} {if (FNR > 1) print}' *FIELDS-simple.txt >> combined_FIELDS2-simple.txt


## keep_samples.txt
# GW10-9 GW10-9
# GW129-24 GW129-24
# GW159-11 GW159-11
# GW167-6 GW167-6
# GW168-1 GW168-1
# GW169-31 GW169-31
# GW30-14 GW30-14
# GW53-5 GW53-5
# GW64-1 GW64-1
# JR19-C3-p22 JR19-C3-p22

# 25-3 and 26-3 were not included.


plink --bfile WGS_153_samples --keep keep_samples.txt --extract extract_hana_variants.txt --keep-allele-order --make-bed --out WGS_153_samples_hana_12_02_2024
plink --bfile WGS_153_samples_hana_12_02_2024 --recodeA --keep-allele-order --out WGS_153_samples_hana_12_02_2024_recodeA 


## Extract VQSR pass variants
for CHR in {1..22} X Y; do
plink --vcf ../CAB_5428.haplotype.recalibrated.decomposed_chr${CHR}.PASS.vcf.gz \
      --extract range gene_regions.bed \
      --make-bed --out WGS_153_samples_VQSR_chr${CHR} \
      --vcf-half-call m --keep-allele-order --double-id
done
#####################
## Extract Emerin  ##
#####################
# Extract ALL EMERIN
plink --chr X --from-bp 154379295 --make-bed --out EMD_WGS_153_samples --to-bp 154381523 --vcf-half-call m --keep-allele-order --double-id --vcf ../CAB_5428.haplotype.recalibrated.decomposed_chrX.vcf.gz
plink --bfile EMD_WGS_153_samples --recodeA --keep-allele-order --out EMD_WGS_153_samples_recodeA