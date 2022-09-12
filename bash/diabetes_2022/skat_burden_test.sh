## Merge all plink files
# plink --bfile MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr1.PASS.decomposed_geno.0.1_hwe.1e-10 --make-bed --merge-list merge_all.list --keep-allele-order --out MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chrALL.PASS.decomposed_geno.0.1_hwe.1e-10


## Extract gnomAD with maf < 0.01 lines
head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt | sed 's/\t/\n/g' | nl
# head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01.txt

for i in {1..22}; do
CHROM="chr${i}"
echo "Doing ${CHROM}"
head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| awk '{print $178"\t"$179"\t"$181"\t"$182"\t"$38"\t"$39"\t"$44}' > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01_${CHROM}.txt	
awk -F'\t' '{if($38 != "." && $38 < 0.01) print $178"\t"$179"\t"$181"\t"$182"\t"$38"\t"$39"\t"$44}' "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHROM}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt" >> /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01_${CHROM}.txt
done