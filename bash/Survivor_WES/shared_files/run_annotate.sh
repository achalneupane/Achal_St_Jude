#!/usr/bin/bash

######################
## Achal Neupane    ##
## Date: 10/10/2023 ##
######################

## Loop over chr to run Survivor_WES_annotate.sh

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=1; \
export INPUT_VCF="${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
# If annotating with gnomAD exome \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/snpEff_round3/"; \
export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/12_10_2023/clinvar.vcf.gz"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/dbnsfp/new_dbNSFP4.4a_variant.${CHR}.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
        -P "${CHR}_annotate" \
        -J "${CHR}_ann" \
        -q "rhel8_gpu" \
        -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
        -gpu "num=1/host" \
        -n ${THREADS} \
        -R "rusage[mem=33192]" \
        "./Survivor_WES_annotate.sh"; \
done;



# Once complete, extract variants:

## 10. Predicted deleterious missense variants: extracted using 1.P_LP_missense.R

## 11. Pathogenic and likely pathogenic variants were extracted based on Clinvar annotation 
cut -d$'\t' -f210 new_chr*.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp_clinvar_12_10_2023_clinvar_12_10_2023.vcf-FIELDS-simple.txt| sort -V | uniq
awk -F'\t' 'NR==1 || ($210 ~ /Pathogenic|Likely_pathogenic/) && $210 !~ /Conflicting|Benign/' \
new_chr*.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp_clinvar_12_10_2023_clinvar_12_10_2023.vcf-FIELDS-simple.txt > all_new_clinvar_P_LP.txt

## 12. High confidence "HC" loss-of-function (excluding NAGNAG site or non-canonincal splice site) variants were extracted
ls chr*.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv| sort -V| xargs cat| grep -v ^##| awk -F'\t' 'NR==1 || ($81 ~ /HC/)' > Loftee_HC_all_chr.txt
# Also remove NAGNAG/ non Canonical splice sites using R
loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)

##  13. All variants were filtered for MAF 0.01 based on gnomAD (all samples) and WES samples. 