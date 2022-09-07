#!/usr/bin/bash


################
## Analysis 2 ##
################
# Yadav:  I think there are some SNPs which are probably <1% in controls but higher in cases and this the overall freq is higher than 1%. Can you remove variants from the summary statistics file that have <1% in controls only and then condition for chr16:25611595:C:T only?

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2
ln -s ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz .
PHENO=../pheno/sjlife_eur_dox_only_pcs.pheno

# extract CMP samples from $PHENO and variants from ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb

awk '{print $1"\t"$2}' $PHENO > samples.list

module load plink/1.90b
plink --chr 16 --from-bp 25361595 --make-bed --out sjlife_CMP --to-bp 25861595 --vcf MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz


awk '{print $2}' ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb > samplesnp.list
# keep samples and variants
ln -s ../finemap/sjlife_CMP.* .
plink --bfile sjlife_CMP --extract samplesnp.list --keep samples.list --keep-allele-order --make-bed --out samplesnp.dat


# plink --bfile sjlife_CMP --freq --keep-allele-order --out sjlife.freq.out
# awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' sjlife.freq.out.frq > sjlife.freq.out.frq_edited1

plink --bfile samplesnp.dat --freq --keep-allele-order --out samplesnp.dat.freq.out
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' samplesnp.dat.freq.out.frq > samplesnp.dat.freq.out.frq_edited1


# Yadav:  I think there are some SNPs which are probably <1% in controls but higher in cases and this the overall freq is higher than 1%. Can you remove variants from the summary statistics file that have <1% in controls only and then condition for chr16:25611595:C:T only?
plink --bfile  samplesnp.dat --keep controls.list --keep-allele-order --make-bed --out samplesnp.dat.controls
plink --bfile samplesnp.dat.controls --freq --keep-allele-order --out samplesnp.dat.controls.freq.out
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' samplesnp.dat.controls.freq.out.frq > samplesnp.dat.controls.freq.out.frq_edited1


plink --bfile sjlife_CMP --extract final.snp.list --keep samples.list --keep-allele-order --make-bed --out samplesnp.dat.final

## Create LD matrix
plink --r2 --bfile samplesnp.dat.final --matrix --out samplesnp_ld_matrix


# # chmod 755 *
# module load dos2unix
# # dos2unix samplesnp.z

##################################
## Run Finemap genomic region 1 ##
##################################
/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 --sss --log --in-files samplesnp --dataset 1

######################
## Run GEC analysis ##
######################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2
module load java
java -jar -Xmx1g /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/gec/gec.jar --var-gec --effect-number --plink-binary  samplesnp.dat.final --genome --out samplesnp_chr16_vars.dat_GEC_test1
# java -jar -Xmx1g gec.jar --var-gec --effect-number --plink-binary samplesnp.dat.final --maf 0.01 --genome --out samplesnp_chr16_vars.dat_GEC_test1ls
# calculate Bonferroni corrected P-value based on Effective_Number in samplesnp_chr16_vars.dat_GEC_test1.sum
Bonferroni_P_chr16 = 0.0002899896 = 0.05/172.42 

###########################
## run GCTA analysis SNP ##
###########################
# module load gcta
# gcta64  --bfile samplesnp.dat.final --chr 16 --cojo-file samplesnp_vars.ma --cojo-slct --cojo-p 0.0002899896 --cojo-wind 100000 --out chr16_cojo
# # gcta64  --bfile samplesnp.dat.final --chr 16 --maf 0.01 --cojo-file samplesnp_TITN_gt_MAF_1_perc_vars.ma --cojo-slct --cojo-p 0.00017925 --cojo-wind 100000 --out TITN_cojo
# # Finally, 4 associated SNPs are selected.

# # Running cojo with selected SNP list 
# gcta64  --bfile samplesnp.dat.final --chr 16 --cojo-file samplesnp_vars.ma --cojo-p 0.0002899896 --cojo-cond chr16_snplist --out chr16_snplist_coho.out
# # gcta64  --bfile samplesnp.dat.final --chr 16 --maf 0.01 --cojo-file samplesnp_vars.ma --cojo-p 0.0002899896 --cojo-cond chr16_snplist --out chr16_25611595_snplist.out

# # # <chr16_snplist>
# # chr16:25554098:G:A
# # chr16:25573423:G:GTA
# # chr16:25609966:C:T
# # chr16:25611595:C:T

# Yadav : "Can you change cojo-p to 5e-08?"
gcta64  --bfile samplesnp.dat.final --chr 16 --cojo-file samplesnp_vars.ma --cojo-slct --cojo-p 5e-08 --cojo-wind 100000 --out chr16_cojo_5e_08
# Finally, 2 associated SNPs are selected...

# chr16:25573423:G:GTA
# chr16:25609966:C:T


# Yadav: "Can you just condition on chr16:25611595:C:T and see if there are additional SNPs with P<5e-08?"
gcta64  --bfile samplesnp.dat.final --chr 16 --cojo-file samplesnp_vars.ma --cojo-p 5e-08 --cojo-cond chr_snplist_chr16_25611595 --out chr16_snplist_chr16_25611595_coho.out




