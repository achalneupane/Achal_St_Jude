# Annotating VCF file with SnpEff. We can add Annovar annotation to the same files or add Annovar annotation separately.

##########################
## 1. SNPEFF Annotation ##
##########################
# download SnpEff and database from: https://pcingola.github.io/SnpEff/download/

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/"; \
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
	"./entrypoint_SnpEff_annotation.sh"; \
done; 



###########################
## 2. Annovar annotation ##
###########################
# # First check avblist hg38_avdblist.txt for the latest releases
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar
# perl annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg38 /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/annovar/annovar/humandb/

# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar 1000g2015aug humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03nontcga humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar esp6500siv2_all humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_exome humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar revel humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar intervar_20180118 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbscsnv11 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic70 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar nci60 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20220320 humandb/ 
# ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar MetaSVM humandb/  # Did not download

# https://annovar.openbioinformatics.org/en/latest/user-guide/startup/
# NOTE: the -operation argument tells ANNOVAR which operations to use for each of the protocols: g means gene-based, gx means gene-based with cross-reference annotation (from -xref argument), r means region-based and f means filter-based.
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar
# Create symlinks to the VCF files on working directory
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf .



for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
bsub \
-P "${CHR}_ANNOVAR" \
-J "${CHR}_ANNOVAR" \
-o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_ANNOVAR.%J" \
-n ${THREADS} \
-R "rusage[mem=8192]" \
"./entrypoint_ANNOVAR_annotation.sh"; \
done;