#!/usr/bin/bash
#################################################
## Helper script to Annotate VCF using ANNOVAR ##
#################################################
## Achal Neupane    ##
## Date: 05/11/2022 ##
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

#################################################################################@@@#################################################
## 1. Create annotation files, we do not need to keep samples. The processed file will therefore be substantially smaller in size. ##
####################################################################################@@@##############################################
CURR_STEP="VCF annotation for ${CHR}"
start=$(${DATE}); echo "[$(display_date ${start})] Starting ${CURR_STEP}"

VCF="${INPUT_VCF}"
MAX_HEADER_LINES=5000
ANNOT_SOURCE="${VCF}"; ANNOT_PROJECT="${VCF%.*}-annot"

# Note: I was getting java version error when I did not reload the modules that I loaded above. So, had to reload them again for each step below. Also, uncomment the commented lines below (steps, all of step 1, 2, 3 and 4) if you are annotating only with ANNOVAR. I have commented these lines because I was addding ANNOVAR annotation to the files that were already annotated by SnpEff and these annotations were already included there. So, these are redundant steps if you already annotated the VCFs with SnpEff. 


#zcat ${ANNOT_SOURCE} | head -${MAX_HEADER_LINES} | grep "^##" > ${ANNOT_PROJECT}.vcf
#zcat ${ANNOT_SOURCE}| grep -v "^##" | cut -f1-8 | awk -F'\t' '!_[$3]++' >> ${ANNOT_PROJECT}.vcf
#sed -i 's/\t\*\t/\t<*:DEL>\t/g' ${ANNOT_PROJECT}.vcf
#echo "DONE trimming the VCF for ${CHR}" >> annotation_step.txt

####################################
## 2. Start annotation with dbSNP ##
####################################
#module load gatk/3.7
#${JAVA} ${JAVAOPTS} -jar ${GATK} \
#   -R ${REF} \
#   -T VariantAnnotator \
#   -V ${ANNOT_PROJECT}.vcf \
#   -L ${ANNOT_PROJECT}.vcf \
#   -o ${ANNOT_PROJECT}.dbSNP155.vcf \
#   --dbsnp /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz
#echo "DONE GATK Annotation with dbSNP for ${CHR}" >> annotation_step.txt

######################################
## 3. Start annotation with clinvar ##
######################################
#module load java/13.0.1
#${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOT_PROJECT}.dbSNP155.vcf > ${ANNOT_PROJECT}.dbSNP155.vcf-clinvar.GRCh38.vcf
#echo "DONE SNPSIFT Annotation with clinvar for ${CHR}" >> annotation_step.txt


#######################
## 4. Adding EXaC db ##
#######################
#${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} ${ANNOT_PROJECT}.dbSNP155.vcf-clinvar.GRCh38.vcf > ${ANNOT_PROJECT}.dbSNP155.vcf-clinvar.ExAC.0.3.GRCh38.vcf
#echo "DONE SNPSIFT Annotation with EXAC for ${CHR}" >> annotation_step.txt


#######################
## 5. Adding ANNOVAR ##
#######################
## Annovar annotation; protocols and operations witll have same number of column tags.
# ANNOTATED="${ANNOT_PROJECT}.dbSNP155.vcf-clinvar.ExAC.0.3.GRCh38.vcf"
ANNOTATED="${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf"
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/table_annovar.pl ${ANNOTATED} \
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/humandb \
-buildver hg38 \
-out ANNOVAR_${ANNOTATED} -remove \
-protocol refGene,1000g2015aug_all,exac03,exac03nontcga,esp6500siv2_all,gnomad_exome,gnomad_genome,dbnsfp42c,intervar_20180118,dbscsnv11,cosmic70,nci60,clinvar_20220320,MetaSVM,revel \
-operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
--nastring . \
-vcfinput
echo "DONE ANNOVAR Annotation for ${CHR}" >> annotation_step.txt

# exitcode=$?
# end=$(${DATE}); echo "[$(display_date ${end})] ${CURR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
# exit ${exitcode}
