## Search variants
# module load parallel
# for chr in {1..22}; do
# export CHR=$(echo $chr)
# echo "Doing CHR${CHR}"
# cat variant_list_uniq.txt | parallel -j30 'cat /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt| grep -w {}' >> matched_positions_snpeff.txt
# done


## Make biallelic
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa WGS_193.GATK4.vcf.gz -Oz -o WGS_193.GATK4.vcf.gz_biallelic.vcf.gz
bcftools index -f -t --threads 4 WGS_193.GATK4.vcf.gz_biallelic.vcf.gz


## Get genotype:
module load plink/1.90b
plink --vcf WGS_193.GATK4.vcf.gz_biallelic.vcf.gz --double-id --vcf-half-call m --keep-allele-order --make-bed --out WGS_193_PLINK


# Missing rsIDs:
# rs104642
# rs8187694


awk '{print $1, $2, $3}' hgTables.txt| tail -n +2|sed 's/ \+/\t/g' > variants_in_vcf.bed
bcftools view -Oz WGS_193.GATK4.vcf.gz_biallelic.vcf.gz --threads 4 -R variants_in_vcf.bed > var_of_interest.vcf.gz


module load tabix
tabix -p vcf var_of_interest.vcf.gz




 ## now anotate only these
 ## SNP eff
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/annotation

unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="var_of_interest.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/annotation/"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
	-P "WGS_193_annotate" \
	-J "WGS_193_annotate" \
	-o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
	-n ${THREADS} \
	-R "rusage[mem=4192]" \
	"./entrypoint_VCFannotation.sh"; 


awk '{split($3,a,";"); for(i in a) if(a[i] ~ /^rs/) $3=a[i]; print}' var_of_interest.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt >  var_of_interest.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple_edited.txt


## add Annovar
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="var_of_interest.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/annotation/"; \
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
"./entrypoint_ANNOVAR_annotation.sh";




module load plink/1.90b
plink --vcf var_of_interest.vcf.gz --double-id --vcf-half-call m --keep-allele-order --make-bed --out var_of_interest_plink
# plink --bfile ${BFILE} --keep-allele-order --allow-no-sex --recode vcf-iid --output-chr MT --out ${VCF}


Work with Rcode: Northwestern_DOC_KO.r; part 2 



## Now recode plink
plink --bfile var_of_interest_plink --recodeA --out var_of_interest_plink_recodeA


Work with Rcode: Northwestern_DOC_KO.r; part 3


###############
## PreQC VCF ##
###############
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/Merge_Intervals_ToChromVCF
bcftools concat WGS_193.haplotype.chr{1..22}.vcf.gz -Oz -o  ./annotation/WGS_193.chr.ALL.tmp.vcf.gz
bcftools view -s WGS_193 WGS_193.chr.ALL.tmp.vcf.gz -Oz -o WGS_193.chr.ALL.vcf.gz
# module load parallel
# for chr in {1..22}; do
# export CHR=$(echo $chr)
# echo "Doing CHR${CHR}"
# cat variant_list_uniq.txt | parallel -j30 'cat /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt| grep -w {}' >> matched_positions_snpeff.txt
# done

## Make biallelic
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa WGS_193.chr.ALL.vcf.gz -Oz -o WGS_193.GATK4.vcf.gz_biallelic.vcf.gz
bcftools index -f -t --threads 4 WGS_193.GATK4.vcf.gz_biallelic.vcf.gz


## Get genotype:
module load plink/1.90b
plink --vcf WGS_193.GATK4.vcf.gz_biallelic.vcf.gz --double-id --vcf-half-call m --keep-allele-order --make-bed --out WGS_193_PLINK


## Now extract variants available in the VCF
# # Variants not found
# rs104642
# rs8187694
awk '{print $1, $2, $3}' hgTables.txt| tail -n +2|sed 's/ \+/\t/g' > variants_in_vcf.bed
bcftools view -Oz WGS_193.GATK4.vcf.gz_biallelic.vcf.gz --threads 4 -R variants_in_vcf.bed > var_of_interest.vcf.gz


module load tabix
tabix -p vcf var_of_interest.vcf.gz




 ## now anotate only these
 ## SNP eff
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/annotation

unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="var_of_interest.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/Merge_Intervals_ToChromVCF/annotation/"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
bsub \
	-P "WGS_193_annotate" \
	-J "WGS_193_annotate" \
	-o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
	-n ${THREADS} \
	-R "rusage[mem=4192]" \
	"./entrypoint_VCFannotation.sh"; 


awk '{split($3,a,";"); for(i in a) if(a[i] ~ /^rs/) $3=a[i]; print}' var_of_interest.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt >  var_of_interest.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple_edited.txt


## add Annovar
unset INPUT_VCF; \
export THREADS=4; \
export INPUT_VCF="var_of_interest.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/Merge_Intervals_ToChromVCF/annotation/"; \
export REF="/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
bsub \
-P "${CHR}_ANNOVAR" \
-J "${CHR}_ANNOVAR" \
-o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_ANNOVAR.%J" \
-n ${THREADS} \
-R "rusage[mem=4192]" \
"./entrypoint_ANNOVAR_annotation.sh";




module load plink/1.90b
plink --vcf var_of_interest.vcf.gz --double-id --vcf-half-call m --keep-allele-order --make-bed --out var_of_interest_plink


Work with Rcode: Northwestern_DOC_KO.r; part 2 



## Now recode plink
plink --bfile var_of_interest_plink --recodeA --out var_of_interest_plink_recodeA


Work with Rcode: Northwestern_DOC_KO.r; part 3



chrX.snp.indel.recalibrated.vcf.gz




