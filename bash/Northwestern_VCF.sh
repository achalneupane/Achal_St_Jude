 # pre VQSR VCF
for file in $(ls *.vcf.gz); do echo "Doing $file" >> bcftools_stats.txt; bcftools stats $file >> bcftools_stats.txt; done
bcftools stats WGS_Northwestern.haplotype.vcf.gz > bcftools_stats_WGS_Northwestern.haplotype.txt

## Post-VQSR
for file in $(ls *.vcf.gz); do echo "Doing $file" >> bcftools_stats_VQSR_PASS.txt; bcftools stats -f PASS $file >> bcftools_stats_VQSR_PASS.txt; done
bcftools stats -f PASS WGS_Northwestern.haplotype.vcf.gz > bcftools_VQSR_PASS_WGS_Northwestern.haplotype.txt


egrep 'Doing|number of SNPs:' bcftools_stats_VQSR_PASS.txt
egrep 'Doing|number of indels:' bcftools_stats_VQSR_PASS.txt


## Get sum of all variants
grep 'number of SNPs:'  bcftools_stats.txt| cut -d$'\t' -f4| paste -sd+ - | bc


 ## Extract variants that PASS VQSR filter
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/WGS_Northwestern_Joint_call_110_samples/WGS_Northwestern/VQSR_ApplyRecalibration



for CHR in {1..22} X Y; do \
    unset VCF; \
    echo "Doing chr${CHR}"; \
    export VCF="chr${CHR}.snp.indel.recalibrated.vcf.gz"; \
    echo -e "**\nDoing ${VCF}"; \
    export MEM=6; \
    export THREADS=4; \
    export INDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/"; \
    export OUT_DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/"; \
    bsub \
    -P "chr${CHR}_edit" \
    -J "chr${CHR}_edit" \
    -o "${OUT_DIR}/logs/${CHR}_s00VCFedit.%J" \
    -n ${THREADS} \
    -R "rusage[mem=6192]" \
    "./Northwestern_edit_VCF.sh"; \
done; 



## Annotate with SNPEff
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/snpEff
ln -s ../../chr*.snp.indel.recalibrated_edited.vcf.gz .
## See VCF_merge_annotate_QC.sh for details about generating reference files and databases

for i in {1..22} X Y; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="${CHR}.snp.indel.recalibrated_edited.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/snpEff/"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
    -P "${CHR}_annotate" \
    -J "${CHR}_annotate" \
    -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
    -n ${THREADS} \
    -R "rusage[mem=8192]" \
    "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/entrypoint_VCFannotation.sh"; \
done; 


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/annovar
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/snpEff/*snp.indel.recalibrated_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf .


## Annotate with Annovar
for i in {1..22} X Y; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="${CHR}.snp.indel.recalibrated_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/annovar/"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
bsub \
-P "${CHR}_ANNOVAR" \
-J "${CHR}_ANNOVAR" \
-o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_ANNOVAR.%J" \
-n ${THREADS} \
-R "rusage[mem=8192]" \
"/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/entrypoint_ANNOVAR_annotation.sh"; \
done;

###########
## Plink ##
###########
## ALL variants
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/plink/ALL

ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/chr*.snp.indel.recalibrated_edited.vcf.gz .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/chr*.snp.indel.recalibrated_edited.vcf.gz.tbi .


for CHR in {1..22} X Y; do \
  unset VCF; unset PLINK_FILE; \
  echo "Doing chr${CHR}"; \
  export VCF="chr${CHR}.snp.indel.recalibrated_edited.vcf.gz"; \
  export PLINK_FILE="chr${CHR}.snp.indel.recalibrated_edited"; \
  export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/plink/ALL/"
  echo -e "**\nConverting ${VCF} \n to Plink \n${PLINK_FILE}\n**"; \
  export MEM=6; \
  export THREADS=4; \
  bsub \
  -P "chr${CHR}_plink" \
  -J "chr${CHR}_plink" \
  -o "${WORKDIR}/logs/${PLINK_FILE}_s01VCFplink.%J" \
  -n ${THREADS} \
  -R "rusage[mem=8192]" \
  "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/plink/entrypoint_VCF_to_Plink_with_basic_QC.sh"; \
done; 


# VQSR PASS variants only
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/plink/VQSR_PASS

ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/chr*.snp.indel.recalibrated_edited.PASS.vcf.gz .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/chr*.snp.indel.recalibrated_edited.PASS.vcf.gz.tbi .


for CHR in {1..22} X Y; do \
  unset VCF; unset PLINK_FILE; \
  echo "Doing chr${CHR}"; \
  export VCF="chr${CHR}.snp.indel.recalibrated_edited.PASS.vcf.gz"; \
  export PLINK_FILE="chr${CHR}.snp.indel.recalibrated_edited.PASS"; \
  export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/plink/VQSR_PASS/"
  echo -e "**\nConverting ${VCF} \n to Plink \n${PLINK_FILE}\n**"; \
  export MEM=6; \
  export THREADS=4; \
  bsub \
  -P "chr${CHR}_plink" \
  -J "chr${CHR}_plink" \
  -o "${WORKDIR}/logs/${PLINK_FILE}_s01VCFplink.%J" \
  -n ${THREADS} \
  -R "rusage[mem=8192]" \
  "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/plink/entrypoint_VCF_to_Plink_with_basic_QC.sh"; \
done; 

## Check MAF 
module load plink/1.90b
for file in $(ls chr*.snp.indel.recalibrated_edited.PASS_geno.0.1_hwe.1e-10.fam); do
bfile=$(basename "$file" .fam)
plink --bfile $bfile --freq --out ${bfile}_freq
done

# compare frequency with /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/1kGP/1000genomes_merged_freq.out.frq


## double check the genotype
bcftools query -s JW10-1 -f '%CHROM %POS  %REF  %ALT [ %GT]\n'  chr1.snp.indel.recalibrated_edited.PASS.vcf.gz > chr1.genotype_JW10-1.txt
bcftools query -s JW6 -f '%CHROM %POS  %REF  %ALT [ %GT]\n'  chr1.snp.indel.recalibrated_edited.PASS.vcf.gz > chr1.genotype_JW6.txt