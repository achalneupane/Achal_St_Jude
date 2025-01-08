## Email on 12/2/2024
Hi Achal,

Happy Monday! Hope you had a great holiday!

I had a question about the WGS for some of the samples we sent you before. Would it be possible to send me the list of variants for the following samples (not sure if 25-3 and 26-3 were included in the ones we sent you)?

19-3
25-3
26-3
GW30
GW64
GW159
GW10
GW53
GW168
GW129
GW167
GW169

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples
module load bcftools

## Fix multiallelic and normalize
bcftools norm -Ou -m -any CAB_5428.haplotype.vcf.gz \
 | bcftools norm -Ou -f /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
 | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' \
 | bcftools view -Oz \
 > CAB_5428.haplotype.recalibrated.vcf.gz


for chr in {1..22} X Y;do
  export chr=${chr}
  bsub -q priority -P HanaWGS -J QC1.chr$chr -eo log/QC1.$chr.err -oo log/QC1.$chr.out -R "rusage[mem=20000]" bash process_WGS_data_perchr.sh $chr
done



#######################
#!/bin/bash

# Processes performed to QC the SJLIFE WGS data provided by Comp. Bio. ["SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714.vcf.gz"]

# Load software packages
module load vcftools/0.1.13
module load bcftools/1.17
#module load vt/2016.11.07
module load plink/1.90b
# module load R/3.4.0
module load R/4.2.2-rhel8
module load bedtools/2.25.0
module load htslib/1.3.1

# Define/create files/directories
INDIR=//research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/

OUTDIR=//research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/

VCFROOT="CAB_5428.haplotype.recalibrated.decomposed"

stats=//research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/stats/

cd $OUTDIR

# Process data per chromosome
chr=$1

echo "Doing chromosome chr${chr}"
# Split data for each chromosome
vcftools --gzvcf ${INDIR}/${VCFROOT}.vcf.gz --chr chr$chr --recode  --stdout | bgzip  > ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz
# bcftools view -S ^${OUTDIR}/chrALL.Survivor_WGS.GATK4180.hg38_missingness.imiss.0.05plus_renamed -r chr$chr -Oz -o ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz ${INDIR}/chrALL.${VCFROOT}.vcf.gz

# Index the VCF file
tabix -pvcf ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz

# Write stats
bcftools stats ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz > ${stats}/${VCFROOT}_chr${chr}.bcftools_stats

# Run basic QC - sequence level
vcftools --gzvcf ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz \
--chr chr$chr \
--keep-filtered PASS \
--minGQ 20 \
--min-meanDP 10 \
--recode \
--stdout \
| bgzip > ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz

# bcftools view -f PASS -i 'MIN(GQ)>20 & MEAN(DP)>10' -Oz -o ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz

# Write stats
bcftools stats ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz > ${stats}/${VCFROOT}_chr${chr}.PASS.bcftools_stats

# Write stats
bcftools stats ${INDIR}/${VCFROOT}_chr${chr}.PASS.decomposed.vcf.gz > ${stats}/${VCFROOT}_chr${chr}.PASS.decomposed.bcftools_stats

tabix -pvcf ${INDIR}/${VCFROOT}_chr${chr}.PASS.decomposed.vcf.gz

# Also write PLINK files for further processing
plink --vcf ${INDIR}/${VCFROOT}_chr${chr}.PASS.decomposed.vcf.gz \
 --keep-allele-order \
 --allow-extra-chr 0 \
 --hwe 1e-10 \
 --make-bed --out ${VCFROOT}_chr${chr}.PASS.decomposed.qced

# Perform LD pruning to obtain independent set of markers
plink --nonfounders \
 --bfile ${VCFROOT}_chr${chr}.PASS.decomposed.qced \
 --exclude range /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/high-LD-regions-GRCh38.txt \
 --geno 0.01 \
 --hwe 0.0001 \
 --maf 0.1 \
 --indep-pairwise 100 25 0.2 \
 --out ${VCFROOT}_chr${chr}.PASS.decomposed.qced_common_pruned
plink --nonfounders \
 --bfile ${VCFROOT}_chr${chr}.PASS.decomposed.qced \
 --extract ${VCFROOT}_chr${chr}.PASS.decomposed.qced_common_pruned.prune.in \
 --make-bed \
 --out ${VCFROOT}_chr${chr}.PASS.decomposed_common_pruned_indep


#######################


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/Hana_request_12_2_2024
module load plink/1.90b

##################
## Extract LMNA ##
##################
# Extract ALL LMNA
plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_WGS_153_samples --to-bp 156140081 --vcf-half-call m --keep-allele-order --vcf CAB_5428.haplotype.recalibrated.decomposed.vcf.gz
plink --bfile LMNA_WGS_153_samples --recodeA --keep-allele-order --out LMNA_WGS_153_samples_recodeA 

#####################
## Extract Emerin  ##
#####################
# Extract ALL EMERIN
plink --chr X --from-bp 154379295 --make-bed --out EMD_WGS_153_samples --to-bp 154381523 --vcf-half-call m --keep-allele-order --vcf CAB_5428.haplotype.recalibrated.decomposed.vcf.gz
plink --bfile EMD_WGS_193 --recodeA --keep-allele-order --out EMD_WGS_193_recodeA