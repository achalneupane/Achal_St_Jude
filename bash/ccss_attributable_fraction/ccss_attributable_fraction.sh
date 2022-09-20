## extract variants from VCF and see which ones are present
module load bcftools/1.9
module load plink/1.90b

ln -s ../../*vcf.gz* .

run extract_variants.py to get variant and PRS score

grep -w chr${CHR} PRS_vars.bed | sort -V > PRS_vars_chr${CHR}.bed
sed -i "s/\r//g" PRS_vars_chr${CHR}.bed
bcftools view -Oz CCSS.GATKv3.4.VQSR_chr${CHR}.PASS.decomposed.ccssid.vcf.gz PRS_vars_chr${CHR}.bed > PRS_chr${CHR}_v3.vcf.gz
echo "Done with extraction for chr${CHR}" >> Status_extract_variant.txt
plink --vcf PRS_chr${CHR}_v3.vcf.gz --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out PRS_chr${CHR}_v3
echo "Done with plink for chr${CHR}" >> Status_extract_variant.txt
