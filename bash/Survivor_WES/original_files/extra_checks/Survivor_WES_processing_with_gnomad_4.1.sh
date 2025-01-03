cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES

hpcf_interactive -n 12 -R "rusage[mem=8000]" -q standard
# ## Convert to biallelic
# bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa MERGED_sorted_sjlife_1_2_zhaoming_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_zhaoming_v2.vcf.gz
# bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' MERGED_biallelic_sorted_sjlife_1_2_zhaoming_v2.vcf.gz -Oz -o MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz
# bcftools index -f -t --threads 4 MERGED_biallelic_sorted_sjlife_1_2_zhaoming_ID_edited.vcf.gz


## Incase we need to rename
awk '{cmd="echo "$0" | sed -e '\''s/.*CCSS-//g; s/^0\\+\\([^0]\\)/\\1/g'\''"; cmd | getline result; close(cmd); print $0"\t"result}' Survivor_WES.samplelist.ccss  > ./biallelic/rename_ccss.txt
cat rename_ccss.txt rename_sjlife.txt > sample_mapping.txt

## <biallelic_renaming_split.sh>
for i in {1..22}; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
unset VCF; \
export THREADS=4; \
export VCF="${CHR}.Survivor_WES.GATK4180.hg38.vcf.gz"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/"; \
bsub \
        -P "${CHR}_biallelic" \
        -J "${CHR}_biallelic" \
        -o "${WORKDIR}/biallelic/logs/${VCF%.vcf*}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "./biallelic_renaming_split.sh"; \
done;


# <biallelic_renaming_split.sh>
#!/usr/bin/bash
module load bcftools
# ## Normalize and make biallelic
bcftools norm -m-any --check-ref -w -f /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa "${WORKDIR}/${VCF}" -Oz -o "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_tmp.vcf.gz"
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_tmp.vcf.gz" -Oz -o "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
bcftools index -f -t --threads 4 "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
# # ## Extract three cohorts
# bcftools view -O z -o \
# "${WORKDIR}/biallelic/ccss/$(basename ${VCF} .vcf.gz)_biallelic_ccss.vcf.gz" \
# -S ${WORKDIR}/biallelic/extract_CCSS.samples.txt \
# "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
# bcftools index -f -t --threads 4 "${WORKDIR}/biallelic/ccss/$(basename ${VCF} .vcf.gz)_biallelic_ccss.vcf.gz"

# bcftools view -O z -o \
# "${WORKDIR}/biallelic/sjlife/$(basename ${VCF} .vcf.gz)_biallelic_sjlife.vcf.gz" \
# -S ${WORKDIR}/biallelic/extract_SJLIFE_survivor.txt \
# "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
# bcftools index -f -t --threads 4 "${WORKDIR}/biallelic/sjlife/$(basename ${VCF} .vcf.gz)_biallelic_sjlife.vcf.gz"

# bcftools view -O z -o \
# "${WORKDIR}/biallelic/sjlife_control/$(basename ${VCF} .vcf.gz)_biallelic_sjlife_control.vcf.gz" \
# -S ${WORKDIR}/biallelic/extract_SJLIFE_survivor_control.txt \
# "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz"
# bcftools index -f -t --threads 4 "${WORKDIR}/biallelic/sjlife_control/$(basename ${VCF} .vcf.gz)_biallelic_sjlife_control.vcf.gz"


## # Rename samples (if we need to rename)
## bcftools reheader -s sample_mapping.txt "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic.vcf.gz" -o "${WORKDIR}/biallelic/$(basename ${VCF} .vcf.gz)_biallelic_renamed.vcf.gz"


## Remove LCR region
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4271055/
# wget https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs38.bed.gz
zcat LCR-hs38.bed.gz| awk '{print $0, "LCR" ++count}'|awk '{sub("chr", "", $1); print $0}'|awk 'BEGIN {OFS="\t"} $1 ~ /^[0-9]+$/ {print $0}' > LCR-hs38.bed

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all
for i in {1..22}; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
unset VCF; \
export THREADS=4; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/"; \
bsub \
        -P "${CHR}_plk" \
        -J "${CHR}_plk" \
        -o "${WORKDIR}/biallelic/plink_all/logs/${VCF%.vcf*}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "./plinkQC.sh"; \
done;


ls chr*.Survivor_WES.GATK4180.hg38_biallelic.bim | sed 's/.bim//'| sort -V > merge_list.txt

# <plinkQC.sh>
#!/usr/bin/bash
## Variant level QC
# Filter based on call rate (<90%)
module load plink/1.90b
#plink --vcf "${WORKDIR}/biallelic/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf.gz" --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic"
plink --bfile "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic" --geno 0.10 --keep-allele-order --make-bed --out "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1"
# Filter based on HWE p-value (<1x10^-15)
plink --bfile "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1" --hwe 1e-15 --keep-allele-order --make-bed --out "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15"
# Exclude variants in low-complexity regions (use BED file)
# bedtools subtract -a filtered_hwe_data.bim -b low_complexity_regions.bed > final_filtered_data.bim
# plink --bfile final_filtered_data --make-bed --out final_filtered_data
plink --bfile "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15" --exclude range LCR-hs38.bed --keep-allele-order --make-bed --out "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed"
# Calculate MAC
module load plink/2.0
plink2 --bfile "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed" --mac 1 --keep-allele-order --make-bed --out "${WORKDIR}/biallelic/plink_all/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1"





# <plinkQC.sh>
#!/usr/bin/bash
## Variant level QC
# Filter based on call rate (<90%)
Biallelic: 2011227
90% call rate: 1832828
HWE: 1816479
LCR removed: 1797175
MAC greater or equal 1: 1796854


## sample level QC
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all
ls chr*.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.bim | sed 's/.bim//'| sort -V  > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15
# plink --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 --indep-pairwise 100 25 0.2 --keep-allele-order --maf 0.05 --make-bed --out chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune


plink \
  --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 \
  --indep-pairwise 100 25 0.2 \
  --keep-allele-order \
  --maf 0.05 \
  --out geno.0.1.hwe.1e-15_indep_pairwise

# plink \
#   --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 \
#   --exclude range high-LD-regions-GRCh38.txt \
#   --indep-pairwise 100 25 0.2 \
#   --keep-allele-order \
#   --maf 0.05 \
#   --out geno.0.1.hwe.1e-15_indep_pairwise

plink --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 \
--keep-allele-order \
--extract geno.0.1.hwe.1e-15_indep_pairwise.prune.in \
--make-bed \
--out chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune

plink --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune --het --out chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune_heterozygosity
plink --bfile chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15 --missing --out chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_missing




## Update IDs
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --update-ids ../../rename_samples.txt --make-bed --keep-allele-order --out chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic --update-ids ../../rename_samples.txt --make-bed --keep-allele-order --out chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated



#######################################################
## Extract Survivorr, control and CCSS exp here from ##
#######################################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all
ls chr*.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.bim | sed 's/.bim//'| sort -V  > merge_list2.txt
plink --merge-list merge_list2.txt --keep-allele-order --make-bed --out chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1


## Survivor; ../extract_SJLIFE_survivor_iid_fid.txt
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --keep-allele-order --keep ../extract_SJLIFE_survivor_iid_fid.txt --make-bed --out Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1

## CCSS exp; ../extract_CCSS.samples_iid_fid.txt
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --keep-allele-order --keep ../extract_CCSS.samples_iid_fid.txt --make-bed --out CCSS_exp/chr.ALL.CCSS_exp_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1

## Control; ../extract_SJLIFE_survivor_control_iid_fid.txt
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1 --keep-allele-order --keep ../extract_SJLIFE_survivor_control_iid_fid.txt --make-bed --out Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1


## Install LOFTEE
https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache
https://wiki.stjude.org/display/CAB/Ensembl+-+VEP
https://wiki.stjude.org/pages/viewpage.action?pageId=40110325

## Revel: https://sites.google.com/site/revelgenomics/downloads; details: https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#revel
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv
## for GRCh37:
## tabix -f -s 1 -b 2 -e 2 new_tabbed_revel.tsv.gz
## for GRCh38:
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz






cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff/revel_loftee
ln -s ../new_chr*.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf .


# for i in {1..22}; do \
# export CHR="chr${i}"; \
# echo "splitting $CHR"; \
# export THREADS=4; \
# export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff/revel_loftee/"; \
# bsub \
#         -P "${CHR}_loftee" \
#         -J "${CHR}_loftee" \
#         -o "${WORKDIR}/logs/${CHR}_loftee.%J" \
#         -n ${THREADS} \
#         -R "rusage[mem=8192]" \
#         "./loftee.sh"; \
# done;

# # <loftee.sh>

# #!/usr/bin/bash
# module load vep/v108
# module load samtools

# cd ${WORKDIR}
# ## Revel
# vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
# --offline --everything --force_overwrite --vcf \
# --fasta /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
# -i new_${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf \
# -o ${CHR}.Survivor_WES_VCF_revel.vcf \
# --assembly GRCh38 \
# --dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/VEP_plugins-release-110 \
# --plugin REVEL,/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/new_tabbed_revel_grch38.tsv.gz


# ## Loftee
# vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
# --offline --everything --force_overwrite \
# --assembly GRCh38 \
# --dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ \
# --fasta /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
# -i ${CHR}.Survivor_WES_VCF_revel.vcf \
# --tab -o ${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.vcf.tsv \
# --plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:false


###############
## 1. SnpEFF ##
###############
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpeff_round3_gnomad4.1
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/*_biallelic.vcf.gz .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/*_biallelic.vcf.gz.tbi .

# # latest dbnsfp 
for CHR in {1..22} M X Y; do 
zcat dbNSFP4.4a_variant.chr${CHR}.gz| bgzip -c > new_dbNSFP4.4a_variant.chr${CHR}.gz
tabix -p vcf new_dbNSFP4.4a_variant.chr${CHR}.gz
done



# for i in {1..22}; do \
# export CHR="chr${i}"; \
# echo "Annotating $CHR"; \
# unset INPUT_VCF; \
# export THREADS=4; \
# export INPUT_VCF="${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf.gz"; \
# export JAVA="java"; \
# export JAVAOPTS="-Xms4g -Xmx30g"; \
# export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff/"; \
# export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
# export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
# export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"
# # export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
# export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
# export EXACDB="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/exac0.3/ExAC.0.3.GRCh38.vcf.gz"; \
# export ONELINEPL="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/scripts/vcfEffOnePerLine.pl"; \
# bsub \
#         -P "${CHR}_annotate" \
#         -J "${CHR}_ann" \
#         -o "${WORKDIR}/logs/${INPUT_VCF%.vcf*}_dbSNP_annotated.%J" \
#         -n ${THREADS} \
#         -R "rusage[mem=8192]" \
#         "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/entrypoint_snpEff_annotation.sh"; \
# done;

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/gnomAD/
## Add gnomad annotation
# wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
# wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/gnomAD/exomes

## Exomes 4.0
for chr in {1..22}; do
echo "Downloading chr ${chr}"
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr${chr}.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr${chr}.vcf.bgz.tbi
done

## Exomes 4.1
for chr in {1..22}; do
echo "Downloading chr ${chr}"
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz.tbi
done



## round 3 gnomad4.1

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
unset INPUT_VCF; \
export THREADS=1; \
export INPUT_VCF="${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf.gz"; \
export JAVA="java"; \
export JAVAOPTS="-Xms4g -Xmx30g"; \
# If annotating with gnomAD exome \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpeff_round3_gnomad4.1/"; \
export REF="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; \
export SNPEFF="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/snpEff.jar"; \
export SNPSIFT="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/SnpSift.jar"; \
# export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/clinvar.vcf.gz"; \
export CLINVAR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/clinvar/12_10_2023/clinvar.vcf.gz"; \
# export GATK="/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar"; \
# export DBNSFPfile="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbNSFP/GRCh38/dbNSFP4.1a.txt.gz"; \
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
        "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/entrypoint_snpEff_annotation_round3_with_gnomad4.1.sh"; \
done;


## to annotate with gnomAD genome
# entrypoint_snpEff_annotation_round3_with_gnomad_gnome.sh


##################################################################
## Helper script to Annotate VCF using snpeff and snpsift tools ##
##################################################################
######################
## Achal Neupane    ##
## Date: 08/16/2024 ##
######################
VERSION="2.0"
#!/usr/bin/bash
# module load gatk/3.7
module load gatk/4.1.8.0
module load vt
module load vcftools
module load bcftools
module load tabix
module load zlib/1.2.5
module load java/13.0.1
module load vep/v108
module load samtools

cd ${WORKDIR}

VCF="${INPUT_VCF}"

MAX_HEADER_LINES=5000
ANNOT_SOURCE="new_${VCF}"; ANNOT_PROJECT="new_${VCF%.*}-annot"

## Adding dbSNP
gatk VariantAnnotator \
   -R ${REF} \
   -V ${VCF} \
   -L ${VCF} \
   -D /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz \
   -O new_${VCF}

echo "DONE GATK Annotation with dbSNP for ${CHR}" >> annotation_step.txt

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

# Adding EXaC db
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v ${EXACDB} ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with EXAC for ${CHR}" >> annotation_step.txt
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf

## Adding clinvar CLNSIG
module load java/13.0.1
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOT_PROJECT}-snpeff-ExAC.0.3.GRCh38.vcf > ${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.vcf
# rm ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf
echo "DONE SNPSIFT Annotation with clinvar for ${CHR}" >> annotation_step.txt


# rm ${ANNOT_PROJECT}-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf
ANNOTATED="${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.vcf"

## Add gnomad
# ## GnomAD with Genome
# ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/gnomAD/gnomad.genomes.v4.0.sites.${CHR}.vcf.bgz \
# ${ANNOTATED} > "${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf"

# ## GnomAD with Exome
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/gnomAD/exomes_4.1/gnomad.exomes.v4.1.sites.${CHR}.vcf.bgz \
"${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.vcf" > "${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf"

ANNOTATED="${ANNOT_PROJECT}-snpeff-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf"


# Revel
vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
--offline --everything --force_overwrite --vcf \
--fasta /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
-i new_${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-ExAC.0.3-clinvar.GRCh38.with.gnomAD.vcf \
-o ${CHR}.Survivor_WES_VCF_revel.vcf \
--assembly GRCh38 \
--dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/VEP_plugins-release-110 \
--plugin REVEL,/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/revel/new_tabbed_revel_grch38.tsv.gz


## Adding dbNSFP
# ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -f '1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,CADD_phred,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate,clinvar_clnsig' -v ${ANNOT_PROJECT}-snpeff.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} dbnsfp -db ${DBNSFPfile} -n -v ${CHR}.Survivor_WES_VCF_revel.vcf > ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf
# rm ${ANNOT_PROJECT}-snpeff.vcf
echo "DONE SNPSIFT annotation with dbnsfp for ${CHR}" >> annotation_step.txt



ANNOTATED=${ANNOT_PROJECT}-snpeff-dbnsfp.vcf

module load java/13.0.1
${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} annotate -v -info CLNSIG ${CLINVAR} ${ANNOTATED} > "$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf"

# # echo "Creating simplified table"
# module load java/13.0.1
# cat ${ANNOTATED} |${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_CADD_phred" "dbNSFP_1000Gp3_AF"  "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_MetaSVM_score" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_pred" "dbNSFP_clinvar_clnsig" "dbNSFP_MutationAssessor_pred"      "dbNSFP_MutationTaster_pred"    "dbNSFP_Polyphen2_HDIV_pred"    "dbNSFP_Polyphen2_HVAR_pred"    "dbNSFP_SIFT_pred"      "dbNSFP_LRT_pred"       "CLNSIG" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_raw" "AF_popmax"> ${ANNOTATED%.*}-FIELDS-simple.txt


ANNOTATED="new_${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp_clinvar_12_10_2023.vcf"
# grep ^## new_chr18.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp.vcf| grep -i dbnsfp| egrep -vi 'gnomad|_AF|_AC ' ''|grep -o 'dbNSFP_[^,]*'| sort | sed 's/.*/"&"/'| tr '\n'
cat ${ANNOTATED}|${ONELINEPL}| ${JAVA} ${JAVAOPTS} -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "dbNSFP_1000Gp3_AF" "dbNSFP_ExAC_AF" "dbNSFP_ExAC_Adj_AF" "dbNSFP_aaalt" "dbNSFP_aapos" "dbNSFP_aaref" "dbNSFP_Aloft_Confidence" "dbNSFP_Aloft_pred" "dbNSFP_Aloft_prob_Dominant" "dbNSFP_Aloft_prob_Recessive" "dbNSFP_Aloft_prob_Tolerant" "dbNSFP_alt" "dbNSFP_AltaiNeandertal" "dbNSFP_APPRIS" "dbNSFP_BayesDel_addAF_pred" "dbNSFP_BayesDel_addAF_rankscore" "dbNSFP_BayesDel_addAF_score" "dbNSFP_BayesDel_noAF_pred" "dbNSFP_BayesDel_noAF_rankscore" "dbNSFP_BayesDel_noAF_score" "dbNSFP_bStatistic" "dbNSFP_bStatistic_converted_rankscore" "dbNSFP_CADD_phred" "dbNSFP_CADD_phred_hg19" "dbNSFP_CADD_raw" "dbNSFP_CADD_raw_hg19" "dbNSFP_CADD_raw_rankscore" "dbNSFP_CADD_raw_rankscore_hg19" "dbNSFP_cds_strand" "dbNSFP_ChagyrskayaNeandertal" "dbNSFP_chr" "dbNSFP_ClinPred_pred" "dbNSFP_ClinPred_rankscore" "dbNSFP_ClinPred_score" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_id" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" "dbNSFP_clinvar_review" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_var_source" "dbNSFP_codon_degeneracy" "dbNSFP_codonpos" "dbNSFP_DANN_rankscore" "dbNSFP_DANN_score" "dbNSFP_Denisova" "dbNSFP_DEOGEN2_pred" "dbNSFP_DEOGEN2_rankscore" "dbNSFP_DEOGEN2_score" "dbNSFP_Eigen_PC_phred_coding" "dbNSFP_Eigen_PC_raw_coding" "dbNSFP_Eigen_PC_raw_coding_rankscore" "dbNSFP_Eigen_phred_coding" "dbNSFP_Eigen_raw_coding" "dbNSFP_Eigen_raw_coding_rankscore" "dbNSFP_Ensembl_geneid" "dbNSFP_Ensembl_proteinid" "dbNSFP_Ensembl_transcriptid" "dbNSFP_FATHMM_converted_rankscore" "dbNSFP_fathmm_MKL_coding_group" "dbNSFP_fathmm_MKL_coding_pred" "dbNSFP_fathmm_MKL_coding_rankscore" "dbNSFP_fathmm_MKL_coding_score" "dbNSFP_FATHMM_pred" "dbNSFP_FATHMM_score" "dbNSFP_fathmm_XF_coding_pred" "dbNSFP_fathmm_XF_coding_rankscore" "dbNSFP_fathmm_XF_coding_score" "dbNSFP_GENCODE_basic" "dbNSFP_genename" "dbNSFP_GenoCanyon_rankscore" "dbNSFP_GenoCanyon_score" "dbNSFP_GERP___NR" "dbNSFP_GERP___RS" "dbNSFP_GERP___RS_rankscore" "dbNSFP_Geuvadis_eQTL_target_gene" "dbNSFP_GM12878_confidence_value" "dbNSFP_GM12878_fitCons_rankscore" "dbNSFP_GM12878_fitCons_score" "dbNSFP_gMVP_rankscore" "dbNSFP_gMVP_score" "dbNSFP_GTEx_V8_gene" "dbNSFP_GTEx_V8_tissue" "dbNSFP_H1_hESC_confidence_value" "dbNSFP_H1_hESC_fitCons_rankscore" "dbNSFP_H1_hESC_fitCons_score" "dbNSFP_hg18_chr" "dbNSFP_hg18_pos_1_based_" "dbNSFP_hg19_chr" "dbNSFP_hg19_pos_1_based_" "dbNSFP_HGVSc_snpEff" "dbNSFP_HGVSc_VEP" "dbNSFP_HGVSp_snpEff" "dbNSFP_HGVSp_VEP" "dbNSFP_HUVEC_confidence_value" "dbNSFP_HUVEC_fitCons_rankscore" "dbNSFP_HUVEC_fitCons_score" "dbNSFP_integrated_confidence_value" "dbNSFP_integrated_fitCons_rankscore" "dbNSFP_integrated_fitCons_score" "dbNSFP_Interpro_domain" "dbNSFP_LINSIGHT" "dbNSFP_LINSIGHT_rankscore" "dbNSFP_LIST_S2_pred" "dbNSFP_LIST_S2_rankscore" "dbNSFP_LIST_S2_score" "dbNSFP_LRT_converted_rankscore" "dbNSFP_LRT_Omega" "dbNSFP_LRT_pred" "dbNSFP_LRT_score" "dbNSFP_M_CAP_pred" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_score" "dbNSFP_MetaLR_pred" "dbNSFP_MetaLR_rankscore" "dbNSFP_MetaLR_score" "dbNSFP_MetaRNN_pred" "dbNSFP_MetaRNN_rankscore" "dbNSFP_MetaRNN_score" "dbNSFP_MetaSVM_pred" "dbNSFP_MetaSVM_rankscore" "dbNSFP_MetaSVM_score" "dbNSFP_MPC_rankscore" "dbNSFP_MPC_score" "dbNSFP_MutationAssessor_pred" "dbNSFP_MutationAssessor_rankscore" "dbNSFP_MutationAssessor_score" "dbNSFP_MutationTaster_AAE" "dbNSFP_MutationTaster_converted_rankscore" "dbNSFP_MutationTaster_model" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationTaster_score" "dbNSFP_MutPred_AAchange" "dbNSFP_MutPred_protID" "dbNSFP_MutPred_rankscore" "dbNSFP_MutPred_score" "dbNSFP_MutPred_Top5features" "dbNSFP_MVP_rankscore" "dbNSFP_MVP_score" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_phastCons100way_vertebrate_rankscore" "dbNSFP_phastCons17way_primate" "dbNSFP_phastCons17way_primate_rankscore" "dbNSFP_phastCons470way_mammalian" "dbNSFP_phastCons470way_mammalian_rankscore" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_phyloP100way_vertebrate_rankscore" "dbNSFP_phyloP17way_primate" "dbNSFP_phyloP17way_primate_rankscore" "dbNSFP_phyloP470way_mammalian" "dbNSFP_phyloP470way_mammalian_rankscore" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HDIV_rankscore" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Polyphen2_HVAR_rankscore" "dbNSFP_Polyphen2_HVAR_score" "dbNSFP_pos_1_based_" "dbNSFP_PrimateAI_pred" "dbNSFP_PrimateAI_rankscore" "dbNSFP_PrimateAI_score" "dbNSFP_PROVEAN_converted_rankscore" "dbNSFP_PROVEAN_pred" "dbNSFP_PROVEAN_score" "dbNSFP_ref" "dbNSFP_refcodon" "dbNSFP_Reliability_index" "dbNSFP_REVEL_rankscore" "dbNSFP_REVEL_score" "dbNSFP_rs_dbSNP" "dbNSFP_SIFT4G_converted_rankscore" "dbNSFP_SIFT4G_pred" "dbNSFP_SIFT4G_score" "dbNSFP_SIFT_converted_rankscore" "dbNSFP_SIFT_pred" "dbNSFP_SIFT_score" "dbNSFP_SiPhy_29way_logOdds" "dbNSFP_SiPhy_29way_logOdds_rankscore" "dbNSFP_SiPhy_29way_pi" "dbNSFP_TSL" "dbNSFP_Uniprot_entry" "dbNSFP_VARITY_ER_LOO_rankscore" "dbNSFP_VARITY_ER_LOO_score" "dbNSFP_VARITY_ER_rankscore" "dbNSFP_VARITY_ER_score" "dbNSFP_VARITY_R_LOO_rankscore" "dbNSFP_VARITY_R_LOO_score" "dbNSFP_VARITY_R_rankscore" "dbNSFP_VARITY_R_score" "dbNSFP_VEP_canonical" "dbNSFP_VEST4_rankscore" "dbNSFP_VEST4_score" "dbNSFP_VindijiaNeandertal" "CLNSIG" "AF_nfe" "AF_afr" "AF_eas" "AF_sas" "AF_raw" "AF" > ${WORKDIR}/"$(basename ${ANNOTATED} .vcf)_clinvar_12_10_2023.vcf-FIELDS-simple2.txt"



# ## Loftee
# vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
# --offline --everything --force_overwrite \
# --assembly GRCh38 \
# --dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ \
# --fasta /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
# -i ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf \
# --vcf -o ${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.vcf \
# --plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:false

## I was not able to split the INFO field for Loftee, so keeping them in a separate folder
mkdir loftee
## Loftee
vep --cache --dir_cache /hpcf/authorized_apps/rhel7_apps/vep/install/108/vep_data/Cache \
--offline --everything --force_overwrite \
--assembly GRCh38 \
--dir_plugins /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/ \
--fasta /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
-i ${ANNOT_PROJECT}-snpeff-dbnsfp.vcf \
--tab -o ${WORKDIR}/loftee/${CHR}.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv \
--plugin LoF,loftee_path:/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/VEP/loftee-grch38/,human_ancestor_fa:false


echo "DONE for ${CHR}" >> annotation_step.txt


#####################
## extract clinvar ##
#####################
cut -d$'\t' -f210 new_chr*.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp_clinvar_12_10_2023_clinvar_12_10_2023.vcf-FIELDS-simple2.txt| sort -V | uniq
awk -F'\t' 'NR==1 || ($210 ~ /Pathogenic|Likely_pathogenic/) && $210 !~ /Conflicting|Benign/' \
new_chr*.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp_clinvar_12_10_2023_clinvar_12_10_2023.vcf-FIELDS-simple2.txt > all_new_clinvar_P_LP.txt

# (base) [aneupane@splprhpc07 snpEff_round3]$ wc -l all_new_clinvar_P_LP.txt
# 99787 all_new_clinvar_P_LP.txt

# extract unique
awk -F'\t' '!seen[$3]++' all_new_clinvar_P_LP.txt > all_new_clinvar_P_LP_unique.txt

#############
### Loftee ##
#############
# grep -v ^## chr22.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv| head -1| sed 's/^#//' > HEADER

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/loftee
ls chr*.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv| sort -V| xargs cat| grep -v ^##| awk -F'\t' 'NR==1 || ($81 ~ /HC/)' > Loftee_HC_all_chr.txt
# awk -F'\t' '!seen[$1]++' Loftee_HC_all_chr.txt > Loftee_HC_all_chr.txt_unique.txt

## Loftee seem to remove the gnomadAD annotation and substitutes it with its own version of gnomAD. To address this, I've reverted to using the gnomADV4 annotation by executing the necessary steps in R.

cd ..


## Extract raw plink
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/clinvar_version_check/gnomad_4.1_12_10_2023
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated* .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated* .

module load plink/1.90b
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated --extract rare_variants_to_extract.txt --keep-allele-order --make-bed -out all_BCC_rare_variants
plink --bfile all_BCC_rare_variants --recodeA --out all_BCC_rare_variants_recodeA