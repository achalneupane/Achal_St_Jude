# Request from Hana emailed on 12/20/2022

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/HANA_Northwestern_request
##################
## Extract LMNA ##
##################
# Extract ALL LMNA
plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_ALL --to-bp 156140081 --vcf-half-call m --keep sample_list.txt --keep-allele-order --vcf chr1.snp.indel.recalibrated_edited.vcf.gz
plink --bfile LMNA_ALL --recodeA --keep-allele-order --out LMNA_ALL_recodeA 

# Extract LMNA_VQSR PASS
plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_ALL_VQSR_PASS --to-bp 156140081 --vcf-half-call m --keep sample_list.txt --keep-allele-order --vcf chr1.snp.indel.recalibrated_edited.PASS.vcf.gz


#####################
## Extract Emerin  ##
#####################
# Extract ALL LMNA
plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_ALL --to-bp 156140081 --vcf-half-call m --keep sample_list.txt --keep-allele-order --vcf chr1.snp.indel.recalibrated_edited.vcf.gz
plink --bfile LMNA_ALL --recodeA --keep-allele-order --out LMNA_ALL_recodeA

# Extract LMNA_VQSR PASS
plink --chr 1 --from-bp 156082573 --make-bed --out LMNA_ALL_VQSR_PASS --to-bp 156140081 --vcf-half-call m --keep sample_list.txt --keep-allele-order --vcf chr1.snp.indel.recalibrated_edited.PASS.vcf.gz


## Now run R script



