####################
## CCSS_Expansion ##
####################
## extract variants from VCF and see which ones are present
module load bcftools/1.9
module load plink/1.90b

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs
# Make VCF biallelic and also edit SNP IDs and header
cat <<\EoF > edit_vcf_header.sh
#!/usr/bin/bash
module load bcftools/1.10.2

cd "${OUT_DIR}"
bcftools norm -m-any --check-ref -w -f /research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa CCSS.vcf.gz -Oz -o CCSS_exp_biallelic_chrALL_tmp1.vcf.gz
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' CCSS_exp_biallelic_chrALL_tmp1.vcf.gz -Oz -o CCSS_exp_biallelic_chrALL_ID_edited_tmp2.vcf.gz
rm MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic.vcf.gz
bcftools reheader -s COMPBIO_ID_TO_CCSS_ID.txt CCSS_exp_biallelic_chrALL_ID_edited_tmp2.vcf.gz > CCSS_exp_biallelic_chrALL_ID_edited.vcf.gz
rm MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_ID_edited.vcf.gz
bcftools index -f -t --threads 4 CCSS_exp_biallelic_chrALL_ID_edited.vcf.gz
EoF


# extract variants from preQC VCF
ln -s ../../CCSS.vcf.gz* .






# run extract_variants.py to get variant and PRS score
for line in $(cat all_cancer_extract_var.txt); do
VAR="$(echo ${line}| tr -d " \t\n\r" )"
# CHR="$(echo $VAR |awk -F':' '{print $1}')"
echo "Doing ${VAR}"
bcftools view CCSS.vcf.gz ${VAR} > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR}.vcf.gz
plink --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --set-missing-var-ids --make-bed --out /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR} 2>&1 | tee -a extract_plink_all.log
done

cd plink_data
# find which ones are missing
for line in $(cat ../all_cancer_extract_var.txt); do
VAR="$(echo ${line}| tr -d " \t\n\r" )"
if [ ! -f "PRS_${VAR}.bim" ] ; then
echo ${VAR} >> failed_plink_files.txt
fi
done

ls *.bim| sort -V | sed 's/\.bim//g'|sed -n '1d;p' > merge_list.list
plink --bfile PRS_chr1:2125052 --merge-list merge_list.list --keep-allele-order --out merged.dat



chr1:109655507 Not found
chr1:145830798 Not found
chr1:172359627 Not found
chr1:51001424 Not found
chr2:9998855 Not found
chr2:217091173 Not found
chr2:39472369 Not found
chr3:141394017 Not found
chr3:49672479 Not found
chr4:105147856 Not found
chr5:344994 Not found
chr5:44508162 Not found
chr5:44619400 Not found
chr5:53383709 Not found
chr5:56366713
chr5:72669180
chr6:151634779
chr6:151701529
chr6:20537614
chr6:81553832
chr7:140243902
chr7:91829875
chr8:123727673
chr8:17930101
chr9:21964883
chr10:22188847
chr10:38234698
chr10:93532430
chr16:3958541
chr17:45134972
chr17:46206492
chr22:44924073
chr1:214274678
chr6:32542614
chr6:32600362
chr12:94126128
chr20:2241198
chr20:34811234
chr5:33946466
chr6:31356838
chr15:28165345
chr6:32641081






















bcftools view -Oz ../CCSS.GATKv3.4.VQSR_chr2.PASS.decomposed.ccssid.vcf.gz chr5:53383700-53383719| zcat | less -S






















#####################
## CCSS _ Original ##
#####################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38


bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' chr21.dose.vcf.gz -Oz -o CCSS_exp_biallelic_chrALL_ID_edited_tmp2.vcf.gz


module load java/17.0.1
module load picard/2.9.4

java -jar /hpcf/apps/picard/install/2.9.4/picard.jar LiftoverVcf \
     I=/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/chr1.dose.vcf.gz \
     O=/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/ccss_org_lifted_over_GRCh38_chr1.vcf \
     CHAIN=/home/aneupane/liftover/hg19ToHg38.over.chain \
     REJECT=/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/rejected_variants_chr1.vcf \
     R=/research_jude/rgs01_jude/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa


