cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/admixture_linux-1.3.0

## SJLIFE

  --geno 0.05
  --hwe 1e-06
  --keep-allele-order
  --maf 0.05
  --make-bed
  --memory 320000
  --merge-list allfiles.txt
  --out SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05


  --bfile SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05
  --extract SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise.prune.in
  --make-bed
  --out dat1


  --bfile dat1
  --make-bed
  --out dat1_renamed
  --update-name rename_sjlife.txt



## 1KG


  --bfile ALL.chr1_GRCh38.genotypes.20170504_biallelic_uniq
  --make-bed
  --nonfounders
  --out ALL.chr1_GRCh38.genotypes.20170504_biallelic_uniq_chrpos
  --update-name ALL.chr1_GRCh38.genotypes.20170504_biallelic_uniq.bim_update_snps



  --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/1kGP/plink/ALL.chr1_GRCh38.genotypes.20170504_biallelic_uniq_chrpos
  --extract to_extract_1kgp.txt
  --make-bed
  --out dat2_1




  --bfile dat2_1
  --make-bed
  --merge-list merge_1kgp.txt
  --out dat2_1kgp





## Merge


  --bfile dat2_1kgp
  --bmerge dat1_renamed
  --make-bed
  --out merged


 --bfile dat2_1kgp_flipped
  --bmerge dat1_renamed
  --make-bed
  --out merged_2


  --bfile dat2_1kgp_flipped_dropped
  --bmerge dat1_renamed
  --make-bed
  --out merged_3


  --bfile dat2_1kgp_flipped_dropped
  --bmerge dat1_renamed_dropped
  --make-bed
  --out merged_4


  --bfile merged_4
  --geno 0.05
  --hwe 1e-06
  --make-bed
  --out final


for K in 1 2 3 4 5; do ../admixture --cv final.bed $K -j8 | tee final_admixture_k${K}.out; done