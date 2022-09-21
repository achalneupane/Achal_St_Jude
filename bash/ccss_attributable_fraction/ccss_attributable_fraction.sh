####################
## CCSS_Expansion ##
####################
## extract variants from VCF and see which ones are present
module load bcftools/1.9
module load plink/1.90b

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs
ln -s ../../*vcf.gz* .

# run extract_variants.py to get variant and PRS score

for line in $(cat all_cancer_extract_var.txt); do
VAR="$(echo ${line}| tr -d " \t\n\r" )"
CHR="$(echo $VAR |awk -F':' '{print $1}')"
echo "Doing ${VAR}"
bcftools view CCSS.GATKv3.4.VQSR_${CHR}.PASS.decomposed.ccssid.vcf.gz ${VAR} > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR}.vcf.gz
plink --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --make-bed --out /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/PRS_${VAR} 2>&1 | tee -a extract_plink_all.log
done

cd plink_data
# find which ones are missing
for line in $(cat ../all_cancer_extract_var.txt); do
VAR="$(echo ${line}| tr -d " \t\n\r" )"
if [ ! -f "PRS_${VAR}.bim" ] ; then
echo ${VAR} >> failed_plink_files.txt
fi
done



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




plink --bfile sjlife_all_PRS_all --bmerge sjlife_all_PRS3 --keep-allele-order --out sjlife_all_PRS_all_final


















#####################
## CCSS _ Original ##
#####################