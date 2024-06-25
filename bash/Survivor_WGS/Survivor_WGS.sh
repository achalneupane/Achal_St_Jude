cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/plink/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4271055/

ln -s ../chr*.Survivor_WGS.GATK4180.hg38.vcf.gz* .

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/plink/
for i in {1..22}; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
unset VCF; \
export VCF="${CHR}.Survivor_WGS.GATK4180.hg38.vcf.gz"
export THREADS=20; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/plink/"; \
bsub \
        -P "${CHR}_plk" \
        -J "${CHR}_plk" \
        -o "${WORKDIR}/logs/${VCF%.vcf*}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=20192]" \
        "./plink.sh"; \
done;




# <plink.sh>
#!/usr/bin/bash
## Variant level QC
# Filter based on call rate (<90%)
module load plink/1.90b
module load bcftools
module load tabix

bcftools norm -m-any --check-ref -w -f /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa ${VCF} --threads 20 -Oz -o "${VCF%.vcf.gz}"_biallelic.vcf.gz
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' "${VCF%.vcf.gz}"_biallelic.vcf.gz --threads 20 -Oz -o "${VCF%.vcf.gz}"_biallelic_edited.vcf.gz
bcftools index -f -t --threads 20 "${VCF%.vcf.gz}"_biallelic_edited.vcf.gz

plink --vcf ${WORKDIR}/"${VCF%.vcf.gz}"_biallelic_edited.vcf.gz --double-id --vcf-half-call m --keep-allele-order --threads 20 --make-bed --out "${WORKDIR}/${CHR}.Survivor_WGS.GATK4180.hg38.vcf_plink"

