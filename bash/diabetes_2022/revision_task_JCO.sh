## Revision task for cindy et al : JCO

SJLIFE_T2D_GWAS_EUR.pheno
SJLIFE_T2D_GWAS_AFR.pheno


## AFR
awk 'NR>1 {print $1, $4}' TOP.AFR.only.with.P.5e-06.and.results.from.meta.variants.txt > TOP.AFR.only.with.P.5e-06.and.results.from.meta.variants_a1 ## list of variants and A1
awk 'NR>1 {print $1}' TOP.AFR.only.with.P.5e-06.and.results.from.meta.variants.txt > TOP.AFR.only.with.P.5e-06.and.results.from.meta.variants ## list of variants only

## EUR
awk 'NR>1 {print $2, $4}' TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants.txt > TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants_a1 ## list of variants and A1
awk 'NR>1 {print $2}' TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants.txt > TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants ## list of variants only

## Sample list 
awk 'NR==1 {print $1, $2; next} {print "0", $2}' ../../AFR_samples.list > AFR_samples.list
awk 'NR==1 {print $1, $2; next} {print "0", $2}' ../../EUR_samples.list > EUR_samples.list

# awk 'NR==1 {print $1, $2; next} {print $1, $2}' ../../AFR_samples.list > AFR_samples.list
# awk 'NR==1 {print $1, $2; next} {print $1, $2}' ../../EUR_samples.list > EUR_samples.list




## extract all variants in AFR TOP.AFR.only.with.P.5e-06.and.results.from.meta.xlsx
module load plink/1.90b







## Afr subset
for chr in {1..22}; do
echo "Doin Chr $chr"
# plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --threads 4 --keep-allele-order --keep AFR_samples.list --make-bed --out plink_chr${chr}.AFR.only
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.preQC_biallelic_renamed_ID_edited.vcf.gz --threads 4 --keep-allele-order --keep AFR_samples.list --make-bed --out plink_chr${chr}.AFR.only
done


## Eur subset
for chr in {1..22}; do
echo "Doin Chr $chr"
# plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --threads 4 --keep-allele-order --keep EUR_samples.list --make-bed --out plink_chr${chr}.EUR.only
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.preQC_biallelic_renamed_ID_edited.vcf.gz --threads 4 --keep-allele-order --keep EUR_samples.list --make-bed --out plink_chr${chr}.EUR.only
done

#############
## African ##
#############
# for file in plink_chr*.AFR.only.bim; do
#     echo "$(basename "$file" .bim)"
# done | sort -V > AFRmerge_list.txt

# plink --merge-list AFRmerge_list.txt --threads 4 --keep-allele-order --make-bed --out All.chr.plink.AFR.only
# rm plink_chr*.AFR.only*

for chr in {1..22}; do
plink --bfile plink_chr${chr}.AFR.only --freq --a1-allele TOP.AFR.only.with.P.5e-06.and.results.from.meta.variants_a1 --extract TOP.AFR.only.with.P.5e-06.and.results.from.meta.variants  --out AFR.Alele1_freq_output_chr${chr}
done

## Discovery
head -1 AFR.Alele1_freq_output_chr1.frq > HEADER
ls AFR.Alele1_freq_output_chr*.frq| sort -V| xargs cat| grep -v CHR > tmp.txt
cat HEADER tmp.txt > AFR.Alele1_freq_output_All.chr

## Replication
for chr in {1..22}; do
plink --bfile plink_chr${chr}.EUR.only --freq --a1-allele TOP.AFR.only.with.P.5e-06.and.results.from.meta.variants_a1 --extract TOP.AFR.only.with.P.5e-06.and.results.from.meta.variants  --out AFR.Alele1_freq_output_chr${chr}_in_EUR
done

head -1 AFR.Alele1_freq_output_chr1_in_EUR.frq > HEADER
ls AFR.Alele1_freq_output_chr*_in_EUR.frq| sort -V| xargs cat| grep -v CHR > tmp.txt
cat HEADER tmp.txt > AFR.Alele1_freq_output_All_in_replication_EUR.chr


## Clumping
## plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --clump chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA --clump-p1 1e-5 --clump-field P --clump-kb 1000 --out clump_output_chr${chr}
## plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --clump clump_output_chr${chr}.clumped --r2 --out LD_results_chr${chr} 
# plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --clump chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA --clump-p1 1e-5 --clump-field P --clump-kb 1000 --r2 --out clump_output_chr${chr}


# # Index variants (EUR)
# chr19:22614854:A:G
# chr2:2434309:AAAATAAAATAAAATAAAATAAAAAC:A
# chr8:50638825:A:G
# # Index variants (AFR)
# chr5:9843018:T:C

## files generated from R script clumping.R
for chr in {1..22}; do
plink --bfile plink_chr${chr}.AFR.only --extract AFR_5e-06_variants_clumping.txt --make-bed --out AFR_discovery_sig_variants_chr${chr}
plink --bfile AFR_discovery_sig_variants_chr${chr} --clump TOP.AFR.only.with.P.5e-06.and.results.txt --clump-snp-field MarkerName --clump-r2 0.1 --r2 --out AFR_LD_clumped_variants_chr${chr}
done

head -1 AFR_LD_clumped_variants_chr1.ld > HEADER
ls AFR_LD_clumped_variants_chr*.ld| sort -V| xargs cat| grep -v CHR_A > tmp.txt
cat HEADER tmp.txt > AFR_LD_clumped_variants_chr_ALL.ld

##############
## European ##
##############
for chr in {1..22}; do
plink --bfile plink_chr${chr}.EUR.only --freq --a1-allele TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants_a1 --extract TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants  --out EUR.Alele1_freq_output_chr${chr}
done

## Discovery
head -1 EUR.Alele1_freq_output_chr1.frq > HEADER
ls EUR.Alele1_freq_output_chr*.frq| grep -v AFR| sort -V| xargs cat| grep -v CHR > tmp.txt
cat HEADER tmp.txt > EUR.Alele1_freq_output_All.chr

## Replication
for chr in {1..22}; do
plink --bfile plink_chr${chr}.AFR.only --freq --a1-allele TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants_a1 --extract TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants  --out EUR.Alele1_freq_output_chr${chr}_in_AFR
done

head -1 EUR.Alele1_freq_output_chr1_in_AFR.frq > HEADER
ls EUR.Alele1_freq_output_chr*_in_AFR.frq| sort -V| xargs cat| grep -v CHR > tmp.txt
cat HEADER tmp.txt > EUR.Alele1_freq_output_All_in_replication_AFR.chr


## Clumping
## plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --clump chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA --clump-p1 1e-5 --clump-field P --clump-kb 1000 --out clump_output_chr${chr}
## plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --clump clump_output_chr${chr}.clumped --r2 --out LD_results_chr${chr} 
# plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --clump chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA --clump-p1 1e-5 --clump-field P --clump-kb 1000 --r2 --out clump_output_chr${chr}


## files generated from R script clumping.R
for chr in {1..22}; do
plink --bfile plink_chr${chr}.EUR.only --extract EUR_5e-06_variants_clumping.txt --make-bed --out EUR_discovery_sig_variants_chr${chr}
plink --bfile EUR_discovery_sig_variants_chr${chr} --clump TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants.txt --clump-snp-field SNP --clump-r2 0.1 --r2 --out EUR_LD_clumped_variants_chr${chr} 
done

head -1 EUR_LD_clumped_variants_chr1.ld > HEADER
ls EUR_LD_clumped_variants_chr*.ld| sort -V| xargs cat| grep -v CHR_A > tmp.txt
cat HEADER tmp.txt > EUR_LD_clumped_variants_chr_ALL.ld