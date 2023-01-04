#!/usr/bin/bash

###############################
## Helper script to edit VCF ##
###############################
######################
## Achal Neupane    ##
## Date: 01/03/2023 ##
######################

VERSION="1.0"

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

# Load module
module load bcftools/1.14

start=$(${DATE}); echo "[$(display_date ${start})] bcftools starting"

if [ ! -f "${VCF}.tbi" ]; then	
echo "Index file missing, Creating index file for ${VCF}"
bcftools index -t --threads ${THREADS} ${VCF}
fi



## Convert to biallelic
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa ${INDIR}/${VCF} -Oz -o ${OUT_DIR}/$(basename "$VCF" .vcf.gz)_v2.vcf.gz
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' ${OUT_DIR}/$(basename "$VCF" .vcf.gz)_v2.vcf.gz -Oz -o ${OUT_DIR}/$(basename "$VCF" .vcf.gz)_edited.vcf.gz
bcftools index -f -t --threads 4 ${OUT_DIR}/$(basename "$VCF" .vcf.gz)_edited.vcf.gz
bcftools view -f PASS ${OUT_DIR}/$(basename "$VCF" .vcf.gz)_edited.vcf.gz -Oz -o ${OUT_DIR}/$(basename "$VCF" .vcf.gz)_edited.PASS.vcf.gz
bcftools index -f -t --threads 4 ${OUT_DIR}/$(basename "$VCF" .vcf.gz)_edited.PASS.vcf.gz ## Extract VQSR PASS only
rm ${OUT_DIR}/$(basename "$VCF" .vcf.gz)_v2.vcf.gz


exit ${exitcode}