#!/usr/bin/bash

######################
## Achal Neupane    ##
## Date: 10/10/2023 ##
######################
VERSION="2.0"
module load gatk/4.1.8.0
module load vt
module load vcftools
module load bcftools
module load tabix
module load zlib/1.2.5
module load java/13.0.1
module load vep/v108
module load samtools

# Loop over each VCF

cd ${WORKDIR}

VCF="${INPUT_VCF}"

MAX_HEADER_LINES=5000
ANNOT_SOURCE="new_${VCF}"; ANNOT_PROJECT="new_${VCF%.*}-annot"

## 1. The script begins by using GATK's VariantAnnotator to annotate variants with dbSNP (V155) information.
gatk VariantAnnotator \
   -R ${REF} \
   -V ${VCF} \
   -L ${VCF} \
   -D /dbSNP/dbSNP_155.GCF_000001405.39.gz \
   -O new_${VCF}
echo "DONE GATK Annotation with dbSNP for ${CHR}" >> annotation_step.txt

## 2. Trim and keep the first 8 columns of VCF. These are the variant IDs we need to annotate.
zcat ${ANNOT_SOURCE}| grep -v "^##" | cut -f1-8 | awk -F'\t' '!_[$3]++' >> ${ANNOT_PROJECT}.vcf
sed -i 's/\t\*\t/\t<*:DEL>\t/g' ${ANNOT_PROJECT}.vcf
echo "DONE trimming the VCF for ${CHR}" >> annotation_step.txt

## 3. Adding snpEff annotation database GRCh38.105. "SnpEff needs a database to perform genomic annotations. There are pre-built databases for thousands of genomes, so chances are that your organism of choice already has a SnpEff database available."
${JAVA} ${JAVAOPTS} -jar ${SNPEFF} -v GRCh38.105  ${ANNOT_PROJECT}.vcf > ${ANNOT_PROJECT}-snpeff.vcf
mv snpEff_genes.txt ${CHR}_snpEff_genes.txt; mv snpEff_summary.html ${CHR}_snpEff_summary.html
echo "DONE SNPeff for ${CHR}" >> annotation_step.txt

## 4. Adding EXaC db release 0.3
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with EXAC for ${CHR}" >> annotation_step.txt

## 5. Adding clinvar database accessed on 12/10/2023: clinvar.vcf.gz
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf > ${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.vcf
echo "DONE SNPSIFT Annotation with clinvar for ${CHR}" >> annotation_step.txt


ANNOTATED="${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.vcf"
## 6. Adding gnomAD (Exome) version 4.0
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/gnomAD/exomes/gnomad.exomes.v4.0.sites.${CHR}.vcf.bgz \
"${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.vcf" > "${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf"
ANNOTATED="${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf"

## 7. Adding Revel annotation plugins-release-110 
vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
--offline --everything --force_overwrite --vcf \
--fasta /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
-i new_${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf \
-o ${CHR}.Survivor_WES_VCF_revel.vcf \
--assembly GRCh38 \
--dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/VEP_plugins-release-110 \
--plugin REVEL,/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/new_tabbed_revel_grch38.tsv.gz


## 8. Adding dbNSFP version 4.4 dbNSFP4.4a_variant.${CHR}.gz
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -n -v ${CHR}.Survivor_WES_VCF_revel.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
# rm ${ANNOT_PROJECT}-snpeff.vcf
echo "DONE SNPSIFT annotation with dbnsfp for ${CHR}" >> annotation_step.txt


## echo "Creating simplified table" from the annotated file. We can extract the columns we need
ANNOTATED="$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf"
cat ${ANNOTATED}|${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF" "dbNSFP_VindijiaNeandertal" "CLNSIG" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_raw" "AF" > ${WORKDIR}/"$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf-FIELDS-simple2.txt"

## 9. Adding Loftee annotation with GRCH38 using VEP plugin version 108
mkdir loftee
vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
--offline --everything --force_overwrite \
--assembly GRCh38 \
--dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ \
--fasta /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
-i ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf \
--tab -o ${WORKDIR}/loftee/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv \
--plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:false

echo "DONE for ${CHR}" >> annotation_step.txt

## 10. Pathogenic and likely pathogenic variants were extracted based on Clinvar annotation 

## 11. High confidence "HC" loss-of-function (excluding NAGNAG site or non-canonincal splice site) variants were extracted

#  12. All variants were filtered for MAF 0.01 based on gnomAD (all samples) and WES samples. 