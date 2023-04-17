module load plink/1.90b
plink --threads 8 --merge-list allfiles.txt --maf 0.05 --hwe 1e-06 --geno 0.05 --make-bed --out SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05

#SJLIFE data
bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE plink --memory 320000 --keep-allele-order --merge-list allfiles.txt --maf 0.05 --hwe 1e-06 --geno 0.05 --make-bed --out SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05
bsub -R "rusage[mem=120000]" -P plinkmergeSJLIFE plink --bfile SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05 --keep-allele-order --allow-extra-chr --exclude range high-LD-regions-GRCh38.txt --geno 0.01 --hwe 0.0001 --maf 0.05 --make-bed --indep-pairwise 100 25 0.2 --out SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise
#1000Genome Data
bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE plink --memory 320000 --keep-allele-order --allow-extra-chr --geno 0.05 --hwe 1e-06 --maf 0.05 --make-bed --out 1000genomes_merged_maf0.05_hwe1e-06_geno0.05 --vcf merged1000genomes.vcf --vcf-half-call h
bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE plink --bfile 1000genomes_merged_maf0.05_hwe1e-06_geno0.05 --memory 320000 --keep-allele-order --allow-extra-chr --make-bed --indep-pairwise 1500 30 0.3 --out 1000genomes_merged_maf0.05_hwe1e-06_geno0.05_indep_pairwise

#Merge 1000 genomes data with SJLIFE data
bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE plink --bfile SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise --bmerge ../../../1kGP/1000genomes_merged_maf0.05_hwe1e-06_geno0.05 --keep-allele-order --allow-extra-chr --make-bed --memory 320000 --out 1000genomes_SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise
'''
1031903 MB RAM detected; reserving 320000 MB for main workspace.
4481 people loaded from
SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise.fam.
2504 people to be merged from
../../../1kGP/1000genomes_merged_maf0.05_hwe1e-06_geno0.05.fam.
Of these, 2504 are new, while 0 are present in the base dataset.
5185765 markers loaded from
SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise.bim.
5120204 markers to be merged from
../../../1kGP/1000genomes_merged_maf0.05_hwe1e-06_geno0.05.bim.
Of these, 5119558 are new, while 646 are present in the base dataset.
Error: 560 variants with 3+ alleles present.
* If you believe this is due to strand inconsistency, try --flip with
  1000genomes_SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise-merge.missnp.
  (Warning: if this seems to work, strand errors involving SNPs with A/T or C/G
  alleles probably remain in your data.  If LD between nearby SNPs is high,
  --flip-scan should detect them.)
* If you are dealing with genuine multiallelic variants, we recommend exporting
  that subset of the data to VCF (via e.g. '--recode vcf'), merging with
  another tool/script, and then importing the result; PLINK is not yet suited
  to handling them.
'''

bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE ../admixture -j8 --supervised final.bed 5
tail -4481 final.5.Q_header2 > final.5.Q_header2_SJLIFE_only

