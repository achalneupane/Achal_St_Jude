cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES


# ## Convert to biallelic
# bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa MERGED_sorted_sjlife_1_2_zhaoming_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_zhaoming_v2.vcf.gz
# bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' MERGED_biallelic_sorted_sjlife_1_2_zhaoming_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz
# bcftools index -f -t --threads 4 MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz

## <biallelic.sh>
#!/usr/bin/bash
module load bcftools
bcftools norm -m-any --check-ref -w -f /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa "${WORKDIR}/${VCF}" -Oz -o "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_tmp.vcf.gz"
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_tmp.vcf.gz" -Oz -o "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
bcftools index -f -t --threads 4 "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset VCF; \
export THREADS=4; \
export VCF="${CHR}.Survivor_WES.GATK4180.hg38.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/"; \
bsub \
        -P "${CHR}_biallelic" \
        -J "${CHR}_biallelic" \
        -o "${WORKDIR}/biallelic/logs/${VCF%.vcf*}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "./biallelic.sh"; \
done;

###############
## 1. SnpEFF ##
###############
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/*_biallelic.vcf.gz .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/*_biallelic.vcf.gz.tbi .


for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff/"; \
export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
# export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
        -P "${CHR}_annotate" \
        -J "${CHR}_ann" \
        -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/entrypoint_snpEff_annotation.sh"; \
done;

## Add gnomad annotation
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi

## Download GRCh38 genomes
for chr in {1..22}; do
echo "Downloading chr ${chr}"
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr${chr}.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr${chr}.vcf.bgz.tbi
done
##################################################################
## Helper script to Annotate VCF using snpeff and snpsift tools ##
##################################################################
######################
## Achal Neupane    ##
## Date: 10/10/2023 ##
######################
# entrypoint_snpEff_annotation.sh
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

VCF="${INPUT_VCF}"
MAX_HEADER_LINES=5000
ANNOT_SOURCE="${VCF}"; ANNOT_PROJECT="${VCF%.*}-annot"

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

## Adding dbNSFP
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -f '1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,CADD_phred,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate,clinvar_clnsig' -v ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
# rm ${ANNOT_PROJECT}-snpeff.vcf
echo "DONE SNPSIFT annotation with dbnsfp for ${CHR}" >> annotation_step.txt

# Adding EXaC db
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with EXAC for ${CHR}" >> annotation_step.txt
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf

## Adding clinvar CLNSIG
module load java/13.0.1
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with clinvar for ${CHR}" >> annotation_step.txt

## Adding dbSNP; note that GATK will load different version of java so will have to load the module again
module load gatk/3.7
${JAVA} ${JAVAOPTS} -jar ${GATK} \
   -R ${REF} \
   -T VariantAnnotator \
   -V ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf \
   -L ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf \
   -o ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf \
   --dbsnp /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz
echo "DONE GATK Annotation with dbSNP for ${CHR}" >> annotation_step.txt
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf
ANNOTATED="${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf"

# echo "Creating simplified table"
module load java/13.0.1
cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_CADD_phred" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_clinvar_clnsig" "dbNSFP_MutationAssessor_pred"	"dbNSFP_MutationTaster_pred"	"dbNSFP_Polyphen2_HDIV_pred"	"dbNSFP_Polyphen2_HVAR_pred"	"dbNSFP_SIFT_pred"	"dbNSFP_LRT_pred"	"CLNSIG"> ${ANNOTATED%.*}-FIELDS-simple.txt

echo "DONE for ${CHR}" >> annotation_step.txt







## Annovar
# <entrypoint>
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


ANNOTATED="${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf"
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/table_annovar.pl ${ANNOTATED} \
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/humandb \
-buildver hg38 \
-out ANNOVAR_${ANNOTATED} -remove \
-protocol refGene,1000g2015aug_all,exac03,exac03nontcga,esp6500siv2_all,gnomad_exome,gnomad_genome,dbnsfp42c,intervar_20180118,dbscsnv11,cosmic70,nci60,clinvar_20220320,MetaSVM,revel \
-operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
-protocol refGene,1000g2015aug_all,exac03,exac03nontcga,esp6500siv2_all,gnomad_exome,gnomad_genome,dbnsfp42c,intervar_20180118,dbscsnv11,cosmic70,nci60,clinvar_20220320,MetaSVM,revel,polyphen2 \
-operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
--nastring . \
-vcfinput
echo "DONE ANNOVAR Annotation for ${CHR}" >> annotation_step.txt
