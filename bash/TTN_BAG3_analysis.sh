#!/bin/bash

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3
# To calculate maf, use sjlife.fam on /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3. Use pheno/sjlife_ttn_bag3.pheno
awk '{print $1"\t"$2}' pheno/sjlife_ttn_bag3.pheno > samples_for_maf.txt
awk '{print $1"\t"$2}' pheno/ccss_org_eur_cardiotoxic_exposed.pheno > samples_for_maf_ccss_org.txt
awk '{print $1"\t"$2}' pheno/ccss_exp_eur_cardiotoxic_exposed.pheno > samples_for_maf_ccss_exp.txt

# USE R code TITN_BAG3.r
# Chunk 1.

# LD calculation with plink; https://zzz.bwh.harvard.edu/plink/ld.shtml; see, "Alternatively, it is possible to add the --matrix option, which creates a matrix of LD values rather than a list: in this case, all SNP pairs are calculated and reported, even for SNPs on different chromosomes"
# Also: https://www.biostars.org/p/268141/
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64
ln -s ../../samples_for_maf.txt


module load plink/1.90b
plink --bfile sjlife --extract samplesnp.list --keep samples_for_maf.txt --keep-allele-order --make-bed --out samplesnp.dat

plink --bfile samplesnp.dat --keep samples_for_maf.txt --freq --keep-allele-order --out sjlife.freq.out

# First removing trailing spaces in the file, then removing chr
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' sjlife.freq.out.frq | awk -F " " '{gsub(/chr/, "", $2); print}' OFS="\t" > sjlife.freq.out.frq_edited1


plink --r2 --bfile samplesnp.dat --matrix --out samplesnp
## samplesnp.list
# chr2:178802888:A:G
# chr2:178803079:C:T
# chr2:178804123:G:A
# ...
# NOTE: The order changes according to the chr and position order in the SNP list. So, the summary statistics needs to be updated to maintain the order as in plink subset. 

# test1.list
# chr2:178313079:G:A
# chr2:178313658:T:C
# chr2:178313672:C:A
# chr2:178313675:A:G
# chr2:178313779:A:G
# chr2:178313920:A:G
# chr2:178314290:C:T
# chr2:178314441:T:C


# test2.list
# chr2:178313079:G:A
# chr2:178313658:T:C
# chr2:178314441:T:C
# chr2:178313920:A:G
# chr2:178313672:C:A
# chr2:178313675:A:G
# chr2:178313779:A:G
# chr2:178314290:C:T




# ## Master file: <samplesnp>
# z;ld;snp;config;cred;log;k;n_samples
# samplesnp.z;samplesnp.ld;samplesnp.snp;samplesnp.config;samplesnp.cred;samplesnp_finemap.log;samplesnp.k;1645
z;ld;snp;config;cred;log;n_samples
samplesnp.z;samplesnp.ld;samplesnp.snp;samplesnp.config;samplesnp.cred;samplesnp_finemap.log;1645



# I was getting errors with the file format, so I had to run dos2unix
# chmod 777 *
# dos2unix samplesnp.z

./finemap_v1.4.1_x86_64 --sss --log --n-causal-snps 2 --in-files samplesnp



# Additionally, .z file shoul have maf in specific form. I was getting this error: "Error : Expected a floating point value greater than 0.0 and smaller or equal than 0.5 in line 3 column 6 of file 'samplesnp.z'!"