cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/annotation/snpeff
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/chr*.Survivor_WGS.GATK4180.hg38.vcf.gz* .

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=2; \
export INPUT_VCF="${CHR}.Survivor_WGS.GATK4180.hg38.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/annotation/snpeff/"; \
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
 -q "rhel8_standard" \
 -J "${CHR}_ann" \
 -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
 -n "${THREADS}" \
 -R "rusage[mem=40GB]" \
 "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/annotation/snpeff/annotate.sh"; \
done;


annotate.sh
##################################################################
## Helper script to Annotate VCF using snpeff and snpsift tools ##
##################################################################
######################
## Achal Neupane    ##
## Date: 10/10/2023 ##
######################
VERSION="2.0"
#!/usr/bin/bash
# module load gatk/3.7
module load gatk/4.1.8.0
module load vt
module load vcftools
module load bcftools
module load tabix
module load zlib/1.2.5
module load java/13.0.1
module load vep/v108
module load samtools

cd ${WORKDIR}


VCF="${INPUT_VCF}"

# bcftools norm -Ou -m -any ${VCF} \
#  | bcftools norm -Ou -f /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
#  | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' \
#  | bcftools view -Oz \
#  > ${VCF}.decomposed.vcf.gz

tabix -p vcf ${VCF}.decomposed.vcf.gz

MAX_HEADER_LINES=5000
ANNOT_SOURCE="new_${VCF}"; ANNOT_PROJECT="new_${VCF%.*}-annot"

## Adding dbSNP
gatk VariantAnnotator \
   -R ${REF} \
   -V ${VCF}.decomposed.vcf.gz \
   -L ${VCF}.decomposed.vcf.gz \
   -D /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz \
   -O new_${VCF}

echo "DONE GATK Annotation with dbSNP for ${CHR}" >> annotation_step.txt

## Start annotating
# zcat ${ANNOT_SOURCE} | head -${MAX_HEADER_LINES} | grep "^##" > ${ANNOT_PROJECT}.vcf
zcat ${ANNOT_SOURCE}| grep -v "^##" | cut -f1-8 | awk -F'\t' '!_[$3]++' >> ${ANNOT_PROJECT}.vcf
sed -i 's/\t\*\t/\t<*:DEL>\t/g' ${ANNOT_PROJECT}.vcf
echo "DONE trimming the VCF for ${CHR}" >> annotation_step.txt
## Adding snpEff annotation; these are genome annotations from ENSEMBL, created from GRCh38/hg38 reference genome sequence
${JAVA} ${JAVAOPTS} -jar ${SNPEFF} -v GRCh38.105  ${ANNOT_PROJECT}.vcf > ${ANNOT_PROJECT}-snpeff.vcf
mv snpEff_genes.txt ${CHR}_snpEff_genes.txt; mv snpEff_summary.html ${CHR}_snpEff_summary.html
# rm ${ANNOT_PROJECT}.vcf
echo "DONE SNPeff for ${CHR}" >> annotation_step.txt

# Adding EXaC db
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with EXAC for ${CHR}" >> annotation_step.txt
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf


module load java/13.0.1

# rm ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf
ANNOTATED="${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf"

## GnomAD with Exome
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/gnomAD/gnomad.genomes.v4.0.sites.${CHR}.vcf.bgz \
"${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf" > "${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.with.gnomAD.vcf"

ANNOTATED="${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.with.gnomAD.vcf"


# Revel
vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
--offline --everything --force_overwrite --vcf \
--fasta /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
-i new_${CHR}.Survivor_WGS.GATK4180.hg38.vcf-annot-snpeff-ExAC.0.3.GRCh38.with.gnomAD.vcf \
-o ${CHR}.Survivor_WES_VCF_revel.vcf \
--assembly GRCh38 \
--dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/VEP_plugins-release-110 \
--plugin REVEL,/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/new_tabbed_revel_grch38.tsv.gz


## Adding dbNSFP
# ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -f '1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,CADD_phred,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate,clinvar_clnsig' -v ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -n -v ${CHR}.Survivor_WES_VCF_revel.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
# rm ${ANNOT_PROJECT}-snpeff.vcf
echo "DONE SNPSIFT annotation with dbnsfp for ${CHR}" >> annotation_step.txt



ANNOTATED=${ANNOT_PROJECT}-snpeff-dbnsfp.vcf

## Adding clinvar CLNSIG
module load java/13.0.1
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOTATED} > "$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf"

# # echo "Creating simplified table"
# module load java/13.0.1
# cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_CADD_phred" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_clinvar_clnsig" "dbNSFP_MutationAssessor_pred"      "dbNSFP_MutationTaster_pred"    "dbNSFP_Polyphen2_HDIV_pred"    "dbNSFP_Polyphen2_HVAR_pred"    "dbNSFP_SIFT_pred"      "dbNSFP_LRT_pred"       "CLNSIG" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_raw" "AF_popmax"> ${ANNOTATED%.*}-FIELDS-simple.txt


ANNOTATED="$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf"
# grep ^## new_chr18.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp.vcf| grep -i dbnsfp| egrep -vi 'gnomad|_AF|_AC ' ''|grep -o 'dbNSFP_[^,]*'| sort | sed 's/.*/"&"/'| tr '\n'
# cat ${ANNOTATED}|${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF" "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_aaalt" "dbNSFP_aapos" "dbNSFP_aaref" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_raw" "AF" > ${WORKDIR}/"$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf-FIELDS-simple2.txt"
OUTFILE="${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple.txt"
cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "dbNSFP_gnomAD_genomes_AF" "dbNSFP_gnomAD_genomes_NFE_AF" "dbNSFP_gnomAD_genomes_AMI_AF" "dbNSFP_gnomAD_genomes_ASJ_AF" "dbNSFP_gnomAD_genomes_MID_AF" "dbNSFP_gnomAD_genomes_FIN_AF" "dbNSFP_gnomAD_genomes_SAS_AF" "dbNSFP_gnomAD_genomes_EAS_AF" "dbNSFP_gnomAD_genomes_AFR_AF" "dbNSFP_gnomAD_genomes_POPMAX_AF" "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_exomes_NFE_AF" "dbNSFP_gnomAD_exomes_ASJ_AF" "dbNSFP_gnomAD_exomes_FIN_AF" "dbNSFP_gnomAD_exomes_SAS_AF" "dbNSFP_gnomAD_exomes_EAS_AF" "dbNSFP_gnomAD_exomes_AFR_AF" "dbNSFP_gnomAD_exomes_POPMAX_AF" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_fin" "AF_raw" "AF" "AF_grpmax" "AF_joint_nfe" "AF_joint_afr" "AF_joint_eas" "AF_joint_sas" "AF_joint_fin" "AF_joint_raw" "AF_joint" "AF_grpmax_joint" > ${WORKDIR}/${$OUTFILE}


## I was not able to split the INFO field for Loftee, so keeping them in a separate folder
mkdir loftee
## Loftee
vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
--offline --everything --force_overwrite \
--assembly GRCh38 \
--dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ \
--fasta /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
-i ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf \
--tab -o ${WORKDIR}/loftee/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv \
--plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:false


echo "DONE for ${CHR}" >> annotation_step.txt







#########

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/loftee
ls chr*.preQC_biallelic_renamed_ID.loftee.tsv| sort -V| xargs cat| grep -v ^##| awk -F'\t' 'NR==1 || ($81 ~ /HC/)' > Loftee_HC_all_chr.txt
awk 'NR > 1 {print $1}' Loftee_HC_all_chr.txt| sort -V | uniq > loftee_unique_HC_variants.txt
## Extract annotated file lines that match HC
for i in {1..22}; do 
export CHR="chr${i}"; 
echo "Working on $CHR"; 
head -1 ../${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > ${CHR}_loftee_HC_variant_extracted_v4.txt
grep -F -f loftee_unique_HC_variants.txt ../${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt >> ${CHR}_loftee_HC_variant_extracted_v4.txt
done

## Extract missense variants
for i in {1..22}; do 
export CHR="chr${i}"; 
echo "Working on $CHR"; 
head -1 ${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > ${CHR}_missense_variant_v4.txt
grep missense ${CHR}.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt >> ${CHR}_missense_variant_v4.txt
done



## P/LP clinvar
cut -d$'\t' -f210 chr*.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt| sort -V | uniq
awk -F'\t' 'NR==1 || ($210 ~ /Pathogenic|Likely_pathogenic/) && $210 !~ /Conflicting|Benign/' \
chr*.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > all_new_clinvar_P_LP.txt

# (base) [aneupane@splprhpc07 snpEff_round3]$ wc -l all_new_clinvar_P_LP.txt
# 99787 all_new_clinvar_P_LP.txt

# extract unique
awk -F'\t' '!seen[$3]++' all_new_clinvar_P_LP.txt > all_new_clinvar_P_LP_unique.txt
