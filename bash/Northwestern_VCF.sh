 # pre VQSR VCF
for file in $(ls *.vcf.gz); do echo "Doing $file" >> bcftools_stats.txt; bcftools stats $file >> bcftools_stats.txt; done
bcftools stats WGS_Northwestern.haplotype.vcf.gz > bcftools_stats_WGS_Northwestern.haplotype.txt

 # Post-VQSR
for file in $(ls *.vcf.gz); do echo "Doing $file" >> bcftools_stats_VQSR_PASS.txt; bcftools stats -f PASS $file >> bcftools_stats_VQSR_PASS.txt; done
bcftools stats -f PASS WGS_Northwestern.haplotype.vcf.gz > bcftools_VQSR_PASS_WGS_Northwestern.haplotype.txt


 egrep 'Doing|number of SNPs:' bcftools_stats_VQSR_PASS.txt
 egrep 'Doing|number of indels:' bcftools_stats_VQSR_PASS.txt


 ## Get sum of all variants
 grep 'number of SNPs:'  bcftools_stats.txt| cut -d$'\t' -f4| paste -sd+ - | bc