#!/usr/bin/bash
##################################################################
## Helper script to Annotate VCF using snpeff and snpsift tools ##
##################################################################
######################
## Achal Neupane    ##
## Date: 05/05/2022 ##
######################
VERSION="1.0"

module load gatk/3.7
module load vt
module load vcftools
module load bcftools
module load tabix
module load vep/v88
module load zlib/1.2.5
module load java/13.0.1

cd ${WORKDIR}

DATE="/bin/date +%s"
display_date () {
    /bin/date -d @${1} +%Y%m%d_%H%M%S
}

date_diff () {
    earlier=${1}
    later=${2}
    diff=$((${later}-${earlier}))
    if [ ${diff} -gt 86400 ]; then
        date -u -d @$((${diff}-86400)) +"%jd%-Hh%-Mm%-Ss"
    else
        date -u -d @${diff} +"%-Hh%-Mm%-Ss"
    fi
}


trap "echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \" && sleep 10s" SIGTERM

CURR_STEP="VCF annotation for ${CHR}"
start=$(${DATE}); echo "[$(display_date ${start})] Starting ${CURR_STEP}"

VCF="${INPUT_VCF}"
MAX_HEADER_LINES=5000 
ANNOT_SOURCE="${VCF}"; ANNOT_PROJECT="${VCF%.*}-annot"
#####################################
## 1. Start annotating with SnpEff ##
#####################################
# create annotation files by removing samples from VCF. This will substantially reduce the size of the VCF files
zcat ${ANNOT_SOURCE} | head -${MAX_HEADER_LINES} | grep "^##" > ${ANNOT_PROJECT}.vcf
zcat ${ANNOT_SOURCE}| grep -v "^##" | cut -f1-8 | awk -F'\t' '!_[$3]++' >> ${ANNOT_PROJECT}.vcf
sed -i 's/\t\*\t/\t<*:DEL>\t/g' ${ANNOT_PROJECT}.vcf
echo "DONE trimming the VCF for ${CHR}" >> annotation_step.txt
## Adding snpEff annotation; these are genome annotations from ENSEMBL, created from GRCh38/hg38 reference genome sequence
${JAVA} ${JAVAOPTS} -jar ${SNPEFF} -v GRCh38.105  ${ANNOT_PROJECT}.vcf > ${ANNOT_PROJECT}-snpeff.vcf
mv snpEff_genes.txt ${CHR}_snpEff_genes.txt; mv snpEff_summary.html ${CHR}_snpEff_summary.html
rm ${ANNOT_PROJECT}.vcf
echo "DONE SNPeff for ${CHR}" >> annotation_step.txt
######################
## 2. Adding dbNSFP ##
######################
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -f '1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,CADD_phred,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate,clinvar_clnsig' -v ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
rm ${ANNOT_PROJECT}-snpeff.vcf
echo "DONE SNPSIFT annotation with dbnsfp for ${CHR}" >> annotation_step.txt

#######################
## 3. Adding EXaC db ##
#######################
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with EXAC for ${CHR}" >> annotation_step.txt
rm ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf 

##############################
## 4. Adding clinvar CLNSIG ##
##############################
module load java/13.0.1
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf
rm ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with clinvar for ${CHR}" >> annotation_step.txt
#####################
## 5. Adding dbSNP ##
##################### 
# note that GATK will load different version of java so will have to load the module again
module load gatk/3.7
${JAVA} ${JAVAOPTS} -jar ${GATK} \
   -R ${REF} \
   -T VariantAnnotator \
   -V ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf \
   -L ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf \
   -o ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf \
   --dbsnp /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz
echo "DONE GATK Annotation with dbSNP for ${CHR}" >> annotation_step.txt
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf
ANNOTATED="${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf"


# echo "Creating simplified table"
module load java/13.0.1
cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_CADD_phred" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_clinvar_clnsig" "CLNSIG"> ${ANNOTATED%.*}-FIELDS-simple.txt
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CURR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
exit ${exitcode}
