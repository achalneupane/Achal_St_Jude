# Request from Hana emailed on 12/20/2022

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/HANA_Northwestern_request
module load plink/1.90b


 zcat chr8.snp.indel.recalibrated.vcf.gz|head -10000| grep "^#CHROM" | cut -f10-|tr '\t' '\n'| grep 19

## Extract LMNA ##
##################
# Extract ALL LMNA
plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_ALL --to-bp 156140081 --vcf-half-call m --keep sample_list.txt --keep-allele-order --vcf chr1.snp.indel.recalibrated_edited.vcf.gz
plink --bfile LMNA_ALL --recodeA --keep-allele-order --out LMNA_ALL_recodeA 

# Extract LMNA_VQSR PASS
plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_ALL_VQSR_PASS --to-bp 156140081 --vcf-half-call m --keep sample_list.txt --keep-allele-order --vcf chr1.snp.indel.recalibrated_edited.PASS.vcf.gz


plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_ALL_FREQ --to-bp 156140081 --vcf-half-call m --keep-allele-order --vcf chr1.snp.indel.recalibrated_edited.vcf.gz
plink --bfile LMNA_ALL_FREQ --freq --out LMNA_ALL_FREQ_result

head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/annovar/ANNOVAR_chr1.snp.indel.recalibrated_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > LMNA_ANNOVAR.txt
grep LMNA  /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/annovar/ANNOVAR_chr1.snp.indel.recalibrated_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt >> LMNA_ANNOVAR.txt

#####################
## Extract Emerin  ##
#####################
# Extract ALL EMERIN
plink --chr X --from-bp 154379295 --make-bed --out EMD_ALL --to-bp 154381523 --vcf-half-call m --keep sample_list.txt --keep-allele-order --vcf chrX.snp.indel.recalibrated_edited.vcf.gz
plink --bfile EMD_ALL --recodeA --keep-allele-order --out EMD_ALL_recodeA

# Extract LMNA_VQSR PASS
plink --chr X --from-bp 154379295 --make-bed --out EMD_ALL_VQSR_PASS --to-bp 154381523 --vcf-half-call m --keep sample_list.txt --keep-allele-order --vcf chrX.snp.indel.recalibrated_edited.PASS.vcf.gz



plink --chr X --from-bp 154379295 --make-bed --out EMD_ALL_FREQ --to-bp 154381523 --vcf-half-call m --keep-allele-order --vcf chrX.snp.indel.recalibrated_edited.vcf.gz
plink --bfile EMD_ALL_FREQ --freq --out EMD_ALL_FREQ_result


head -1  /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/annovar/ANNOVAR_chrX.snp.indel.recalibrated_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > EMD_ANNOVAR.txt
grep EMD  /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/annotation/annovar/ANNOVAR_chrX.snp.indel.recalibrated_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt >> EMD_ANNOVAR.txt

## Now run R script Hana_12_20_2022_email_request_for_variant_extraction.R





## (9/18/2024)
# Hi Achal,
# Hope all is good.
# I had a question about the WGS we sent before. As far as I remember, we included 19-3 cell line as well. Would it be possible to send me the variants in LMNA gene for that cell line? You previously kindly provided me with the attached excel file. Now we need the exact same thing for 19-3 cell line.
# Thanks a lot,
# Best,
# Hana

cd  /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/

bcftools norm -Ou -m -any WGS_193.GATK4.vcf.gz \
 | bcftools norm -Ou -f /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/hg38.fa \
 | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' \
 | bcftools view -Oz \
 > WGS_193.GATK4.recalibrated.decomposed.vcf.gz


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/HANA_Northwestern_request
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193
## Extract LMNA ##
##################
# Extract ALL LMNA
plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_WGS_193 --to-bp 156140081 --vcf-half-call m --keep-allele-order --vcf WGS_193.GATK4.recalibrated.decomposed.vcf.gz
plink --bfile LMNA_WGS_193 --recodeA --keep-allele-order --out LMNA_WGS_193_recodeA 


#####################
## Extract Emerin  ##
#####################
# Extract ALL EMERIN
plink --chr X --from-bp 154379295 --make-bed --out EMD_WGS_193 --to-bp 154381523 --vcf-half-call m --keep-allele-order --vcf WGS_193.GATK4.recalibrated.decomposed.vcf.gz
plink --bfile EMD_WGS_193 --recodeA --keep-allele-order --out EMD_WGS_193_recodeA

