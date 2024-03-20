## Plink subset
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/final_data_submission
## Eur subset
awk 'NR==1 {print $1, $2; next} {print "0", $2}' EUR_samples.list > EUR_samples.list_edited
for chr in {1..22}; do
echo "Doin Chr $chr"
# plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --threads 4 --keep-allele-order --keep EUR_samples.list --make-bed --out plink_chr${chr}.EUR.only
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.preQC_biallelic_renamed_ID_edited.vcf.gz --threads 4 --keep-allele-order --keep ../EUR_samples.list_edited --make-bed --out SJLIFE_EUR_chr${chr}
done

## Afr subset
awk 'NR==1 {print $1, $2; next} {print "0", $2}' AFR_samples.list > AFR_samples.list_edited
for chr in {1..22}; do
echo "Doin Chr $chr"
# plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${chr}.PASS.decomposed_geno.0.1_hwe.1e-10 --threads 4 --keep-allele-order --keep AFR_samples.list --make-bed --out plink_chr${chr}.AFR.only
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.preQC_biallelic_renamed_ID_edited.vcf.gz --threads 4 --keep-allele-order --keep ../AFR_samples.list_edited --make-bed --out SJLIFE_AFR_chr${chr}
done


## Compress files
tar -czvf plink_data_AFR_compressed.tar.gz plink_data_AFR
tar -czvf plink_data_EUR_compressed.tar.gz plink_data_EUR