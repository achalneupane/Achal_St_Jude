#!/bin/bash

mkdir HDL
cp AFR_HDLsumary
cp EUR_Sumaary


EUR=DYSLPDM_chltotl_rank_chr*.assoc.linear


## Process PLINK result files for meta-analysis using METAL; add A2 and BETA
# SJLIFE
# Find the header file and save it to a variable
header_file=$(ls DYSLPDM_chltotl_rank_chr*.assoc.linear | head -n 1)
# Concatenate all files, skipping the header for each file except the first
cat $header_file > all_chr_EUR  # Copy the header from the first file
cat $(ls DYSLPDM_chltotl_rank_chr*.assoc.linear | sort -V) | grep -v CHR >> all_chr_EUR

tr '\t' ' ' < all_chr_EUR | sed 's/  */ /g' > DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean
rm all_chr_EUR
{ head -n 1 DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean && tail -n +2 DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean | sort -k12,12g; } > DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted


AFR=DYSLPDM_chltotl_rankeur_afr_chr*.assoc.linear

header_file=$(ls DYSLPDM_chltotl_rankeur_afr_chr*.assoc.linear | head -n 1)
# Concatenate all files, skipping the header for each file except the first
cat $header_file > all_chr_AFR  # Copy the header from the first file
cat $(ls DYSLPDM_chltotl_rankeur_afr_chr*.assoc.linear | sort -V) | grep -v CHR >> all_chr_AFR

tr '\t' ' ' < all_chr_AFR | sed 's/  */ /g' > DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean
rm all_chr_AFR
{ head -n 1 DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean && tail -n +2 DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean | sort -k12,12g; } > DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted




awk 'BEGIN {OFS="\t"} NR==1 {print $0, "OR"} NR>1 {print $0, exp($7)}' DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted > tmp
awk 'BEGIN {OFS=" "} {print $1, $2, $3, $4, $5, $6, $13, $14, $8, $9, $10, $11, $12}' tmp > DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted

awk 'BEGIN {OFS="\t"} NR==1 {print $0, "OR"} NR>1 {print $0, exp($7)}' DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted > tmp
awk 'BEGIN {OFS=" "} {print $1, $2, $3, $4, $5, $6, $13, $14, $8, $9, $10, $11, $12}' tmp > DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted
rm tmp

# EUR
awk '{split($2, a, ":"); print $0, a[3], a[4]}' DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted \
| awk '{if($4==$13) A2=$14; else A2=$13; BETA=log($7); print $0, A2, BETA}' \
| awk 'BEGIN{print "CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA"}FNR>1{print}' \
> DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal
# AFR
awk '{split($2, a, ":"); print $0, a[3], a[4]}' DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted \
| awk '{if($4==$13) A2=$14; else A2=$13; BETA=log($7); print $0, A2, BETA}' \
| awk 'BEGIN{print "CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P REF ALT A2 BETA"}FNR>1{print}' \
> DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal



## Identify variants that match between the 2 datasets, based on chr:pos and then based on both A1 and A2
awk 'NR==FNR{a[$1":"$3]=$0;next}($1":"$3 in a){print $0, a[$1":"$3]}' DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal \
DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal| awk '($4==$20 || $4==$31) && ($15==$20 || $15==$31)' > common_variants_among_2_datasets.txt

## Fix space
sed 's/[[:blank:]]\+/ /g' common_variants_among_2_datasets.txt > common_variants_among_2_datasets_updated.txt

# Then prepare final files for each dataset for METAL
awk '!a[$1":"$3]++' common_variants_among_2_datasets_updated.txt | cut -d' ' -f1-16 | awk '{if(FNR>1) $2=$1":"$3; print}' \
> DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal.common
awk '!a[$1":"$3]++' common_variants_among_2_datasets_updated.txt | cut -d' ' -f17-32 | awk '{if(FNR>1) $2=$1":"$3; print}' \
> DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal.common



## Prepare the config file for METAL
# run_metal.script

## Run METAL
~/bin/generic-metal/metal run_metal.script
SCHEME STDERR
STDERR SE

MARKER SNP
ALLELE A1 A2
EFFECT BETA
PVALUE P
PROCESS DYSLPDM_chltotl_rank_chrALL.assoc.linear.clean.Psorted.formetal.common

PROCESS DYSLPDM_chltotl_rankeur_afr_chrALL.assoc.linear.clean.Psorted.formetal.common

OUTFILE Meta_DYSLPDM_chltotl_EUR_AFR_fixed_ .tbl

ANALYZE HETEROGENEITY

QUIT



## A total of 7319447 variants present in both datasets were meta-analyzed; process the results further
module load plink/1.90b
DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq"
# for CHR in {1..22}; do
# plink --bfile MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/SJLIFE_EURsamples.txt --freq --out "${DIR}/dyslypidemia_EUR_frequency_chr${CHR}"
# awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' "${DIR}/dyslypidemia_EUR_frequency_chr${CHR}.frq" | awk -F " " '{gsub(/chr/, "", $2); print}' OFS="\t" > "${DIR}/dyslypidemia_EUR_frequency_chr${CHR}_edited1.frq"
# done

# for CHR in {1..22}; do
# plink --bfile MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/SJLIFE_AFRsamples.txt --freq --out "${DIR}/dyslypidemia_AFR_frequency_chr${CHR}"
# awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' "${DIR}/dyslypidemia_AFR_frequency_chr${CHR}.frq" | awk -F " " '{gsub(/chr/, "", $2); print}' OFS="\t" > "${DIR}/dyslypidemia_AFR_frequency_chr${CHR}_edited1.frq"
# done

# ## AFR
# head -1 dyslypidemia_AFR_frequency_chr1_edited1.frq > HEADER
# ls dyslypidemia_AFR_frequency_chr*_edited1.frq|sort -V| xargs cat | grep -v ^CHR> tmp.frq
# cat HEADER tmp.frq > ALL_CHR_dyslypidemia_AFR_frequency_edited1.frq
# rm tmp.frq 
# rm dyslypidemia_AFR_frequency_chr*_edited1.frq
# rm dyslypidemia_AFR_frequency_chr*


# ## EUR
# head -1 dyslypidemia_EUR_frequency_chr1_edited1.frq > HEADER
# ls dyslypidemia_EUR_frequency_chr*_edited1.frq|sort -V| xargs cat | grep -v ^CHR> tmp.frq
# cat HEADER tmp.frq > ALL_CHR_dyslypidemia_EUR_frequency_edited1.frq
# rm tmp.frq 
# rm dyslypidemia_EUR_frequency_chr*_edited1.frq
# rm dyslypidemia_EUR_frequency_chr*

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq
DIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq"
ln -s ../MTAG/* .

module load plink/1.90b

## We will need to extract variants from one of the studies. Since we will ultimately be using metal to combine all studies, we just need to extract the overlapping variants from the studies. Note, they all have almost same variants, so I am just using this file to extract the variants

# awk 'NR==FNR{a[$1":"$4]=$0;next}($1":"$3 in a){print $0, a[$1":"$3]}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/ALL_bim.txt DYSLPDM_eur_chltot_chrALL.assoc.linear.clean.Psorted.formetal.common  | awk '($4==$21 || $4==$22) && ($15==$21 || $15==$22)'> variants_to_extract

# awk 'NR==FNR{a[$1":"$4]=$0;next}($1":"$3 in a){print $0, a[$1":"$3]}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/ALL_bim.txt DYSLPDM_eur_chltot_chrALL.assoc.linear.clean.Psorted.formetal.common  | awk '(($4==$21 || $4==$22) && ($13==$21 || $13==$22))||(($4==$21 || $4==$22) && ($15==$21 || $15==$22))'> variants_to_extract

awk 'NR==FNR{a[$1":"$3]=$0;next}($1":"$4 in a){print $0, a[$1":"$4]}' DYSLPDM_eur_chltot_chrALL.assoc.linear.clean.Psorted.formetal.common /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/ALL_bim.txt   | awk '(($5==$19 || $5==$20) && ($6==$19 || $6==$20))||(($5==$20 || $5==$21) && ($6==$20 || $6==$21))' > variants_to_extract1

awk '{print $2}' variants_to_extract1 > variants_to_extract2


#########
## EUR ##
#########
for CHR in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/SJLIFE_EURsamples.txt --extract variants_to_extract2 --make-bed --out "${DIR}/dyslypidemia_EUR_chr${CHR}"
awk -F'\t' '{split($2, a, ":"); $2 = a[1] ":" a[2]}1' "${DIR}/dyslypidemia_EUR_chr${CHR}.bim" > "${DIR}/tmp.bim"
mv "${DIR}/tmp.bim" "${DIR}/dyslypidemia_EUR_chr${CHR}.bim"
done

# Generate a list of file names for chromosomes 1 to 22
for i in {1..22}; do echo "dyslypidemia_EUR_chr${i}.bed dyslypidemia_EUR_chr${i}.bim dyslypidemia_EUR_chr${i}.fam"; done > merge_list.txt
# Use Plink to merge the files
plink --merge-list merge_list.txt --make-bed --out merged_dyslipidemia_EUR

TYPES=("DYSLPDM_eur_chltot" "DYSLPDM_eur_hdl" "DYSLPDM_eur_ldl" "DYSLPDM_eur_nonhdl" "DYSLPDM_eur_trigly")
for TYPE in "${TYPES[@]}"; do
    awk 'NR > 1 {print "chr"$2, $4}' "${TYPE}_chrALL.assoc.linear.clean.Psorted.formetal.common" > "allele_${TYPE}.txt"
    plink --bfile "${DIR}/merged_dyslipidemia_EUR" --freq --a1-allele "${DIR}/allele_${TYPE}.txt" --out "${DIR}/${TYPE}_frequency"
    awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' "${DIR}/${TYPE}_frequency.frq" | awk -F " " '{gsub(/chr/, "", $2); print}' OFS="\t" > "${DIR}/${TYPE}_frequency_edited.frq"
done


#########
## AFR ##
#########
for CHR in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed --keep /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/SJLIFE_AFRsamples.txt --extract variants_to_extract2 --make-bed --out "${DIR}/dyslypidemia_AFR_chr${CHR}"
awk -F'\t' '{split($2, a, ":"); $2 = a[1] ":" a[2]}1' "${DIR}/dyslypidemia_AFR_chr${CHR}.bim" > "${DIR}/tmp.bim"
mv "${DIR}/tmp.bim" "${DIR}/dyslypidemia_AFR_chr${CHR}.bim"
done

# Generate a list of file names for chromosomes 1 to 22
for i in {1..22}; do echo "dyslypidemia_AFR_chr${i}.bed dyslypidemia_AFR_chr${i}.bim dyslypidemia_AFR_chr${i}.fam"; done > merge_listAFR.txt
# Use Plink to merge the files
plink --merge-list merge_listAFR.txt --make-bed --out merged_dyslipidemia_AFR

TYPES=("DYSLPDM_afr_chltot" "DYSLPDM_afr_hdl" "DYSLPDM_afr_ldl" "DYSLPDM_afr_nonhdl" "DYSLPDM_afr_trigly")
for TYPE in "${TYPES[@]}"; do
    awk 'NR > 1 {print "chr"$2, $4}' "${TYPE}_chrALL.assoc.linear.clean.Psorted.formetal.common" > "allele_${TYPE}.txt"
    plink --bfile "${DIR}/merged_dyslipidemia_AFR" --freq --a1-allele "${DIR}/allele_${TYPE}.txt" --out "${DIR}/${TYPE}_frequency"
    awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' "${DIR}/${TYPE}_frequency.frq" | awk -F " " '{gsub(/chr/, "", $2); print}' OFS="\t" > "${DIR}/${TYPE}_frequency_edited.frq"
done


conda create --name python2.7 python=2.7
conda env list
conda activate /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/envs/python2.7
conda install scipy.optimize

# # /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/mtag.py --sumstats eur_chltot.shared.txt,eur_ldl.shared.txt,eur_hdl.shared.txt,eur_trigly.shared.txt,eur_nonhdl.shared.txt --snp_name snpid --z_name z --ld_ref_panel /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/ld_ref_panel/eur_w_ld_chr_snp --out ./EUR_chltot_ldl_hdl_trigly_nonhdl_mtag --force --n_min 0.0 --stream_stdout
# # Warning: The mean chi2 statistic of trait 5 is less 1.02 - MTAG estimates may be unstable.
# # Dropped 1022520 SNPs due to strand ambiguity, 5549574 SNPs remain in intersection after merging trait1
# # Dropped 1269 SNPs due to inconsistent allele pairs from phenotype 2. 5547254 SNPs remain.
# # Dropped 0 SNPs due to strand ambiguity, 5547254 SNPs remain in intersection after merging trait2
# # Dropped 1961 SNPs due to inconsistent allele pairs from phenotype 3. 5543833 SNPs remain.
# # Dropped 0 SNPs due to strand ambiguity, 5543833 SNPs remain in intersection after merging trait3
# # Dropped 742 SNPs due to inconsistent allele pairs from phenotype 4. 5542444 SNPs remain.
# # Dropped 0 SNPs due to strand ambiguity, 5542444 SNPs remain in intersection after merging trait4
# # Dropped 110 SNPs due to inconsistent allele pairs from phenotype 5. 5542263 SNPs remain.
# # Dropped 0 SNPs due to strand ambiguity, 5542263 SNPs remain in intersection after merging trait5
# # ... Merge of GWAS summary statistics complete. Number of SNPs:   5542263
# # Using 5542263 SNPs to estimate Omega (0 SNPs excluded due to strand ambiguity)
# # Estimating sigma..
# # After merging with reference panel LD, 0 SNPs remain.
# # Traceback (most recent call last):
# #   File "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/mtag.py", line 1577, in <module>
# #     mtag(args)
# #   File "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/mtag.py", line 1358, in mtag
# #     args.sigma_hat = estimate_sigma(DATA[not_SA], args)
# #   File "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/mtag.py", line 472, in estimate_sigma
# #     rg_results =  sumstats_sig.estimate_rg(args_ldsc_rg, Logger_to_Logging())
# #   File "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/ldsc_mod/ldscore/sumstats.py", line 423, in estimate_rg
# #     M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(args, log, None, alleles=True, dropna=True,sumstats=p1)
# #   File "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/ldsc_mod/ldscore/sumstats.py", line 251, in _read_ld_sumstats
# #     sumstats = _merge_and_log(ref_ld, sumstats, 'reference panel LD', log)
# #   File "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/ldsc_mod/ldscore/sumstats.py", line 236, in _merge_and_log
# #     raise ValueError(msg.format(N=len(sumstats), F=noun))
# # ValueError: After merging with reference panel LD, 0 SNPs remain.
# # Analysis terminated from error at Tue Feb 27 15:00:13 2024
# # Total time elapsed: 6.0m:35.79s

# # # https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# # # https://github.com/JonJala/mtag/wiki/Tutorial-1:-The-Basics

# # --ld_ref_panel FOLDER_PATH
# #                         Specify folder of the ld reference panel (split by
# #                         chromosome) that will be used in the estimation of the
# #                         error VCV (sigma). This option is passed to --ref-ld-
# #                         chr and --w-ld-chr when running LD score regression.
# #                         The default is to use the reference panel of LD scores
# #                         computed from 1000 Genomes European subjects
# #                         (eur_w_ld_chr) that is included with the distribution
# #                         of MTAG


# # uses rsID, so modfying in another copy of ld file
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/ld_ref_panel
# cp -R eur_w_ld_chr eur_w_ld_chr_snp
# for file in *ldscore.gz; do
# 	echo "doing ${file}"
# gunzip -c ${file} | awk 'NR==1 {print $0; next} {$2 = $1":"$3; print}'| gzip > tmp.gz
# mv tmp.gz ${file}
# done



## Now trying again with our own ld panel
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG
 
## EUR
# /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldsc/ldsc.py --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/merged_dyslipidemia_EUR --l2 --ld-wind-cm 1 --yes-really --out merged_dyslipidemia_EUR_LD
# # (python2.7) [aneupane@splprhpc09 MTAG]$ /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldsc/ldsc.py --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/merged_dyslipidemia_EUR --l2 --ld-wind-cm 1 --out merged_dyslipidemia_EUR_LD
module load plink/1.90b
# run per chromosome
# for CHR in {1..22}; do
# plink --bfile  /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/dyslypidemia_EUR_chr${CHR} --cm-map --make-bed --out dyslypidemia_EUR_chr${CHR}_ldscore
# /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldsc/ldsc.py --bfile dyslypidemia_EUR_chr${CHR}_ldscore --l2 --ld-wind-cm 1 --yes-really --out dyslipidemia_EUR_chr${CHR}_LD
# done

https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial

https://alkesgroup.broadinstitute.org/Eagle/
## download genetic map from
wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz

##############
## European ##
##############

# for CHR in {1..22}; do
# # /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldsc/ldsc.py --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq//dyslypidemia_EUR_chr${CHR} --l2 --ld-wind-kb 1000 --yes-really --out $PWD/ldscore/dyslipidemia_EUR_chr${CHR}_LD
# zcat /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/genetic_map_hg38_withX.txt.gz |grep -w ${CHR}| awk -v chr="${CHR}" '$1 == chr'|awk '{print $2, $3, $4}' > $PWD/ldscore/genetic_map_chr${CHR}_combined_hg38.txt
# plink --bfile  /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/dyslypidemia_EUR_chr${CHR} --cm-map $PWD/ldscore/genetic_map_chr@_combined_hg38.txt --make-bed --out $PWD/ldscore/dyslypidemia_EUR_chr${CHR}_ldscore
# /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldsc/ldsc.py --bfile $PWD/ldscore/dyslypidemia_EUR_chr${CHR}_ldscore --l2 --ld-wind-cm 1 --yes-really --out $PWD/ldscore/dyslipidemia_EUR_chr${CHR}_LD
# done


# # OR liftover: https://mathgen.stats.ox.ac.uk/impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html
# wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
# tar -xzvf 1000GP_Phase3.tgz
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3

module load plink/1.90b

for CHR in {1..22}; do
awk -v CHR=${CHR} '{print CHR, $0}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/genetic_map_chr${CHR}_combined_b37.txt > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/tmpfile

awk 'BEGIN {OFS="\t"} NR>1 {print "chr"$1, $2, $2, $3, $4, $5}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/tmpfile > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/modified_genetic_map_chr${CHR}_combined_b37.bed

/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/liftover/liftOver /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/modified_genetic_map_chr${CHR}_combined_b37.bed /home/aneupane/liftover/hg19ToHg38.over.chain /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/genetic_map_chr${CHR}_combined_hg38_mapped.bed /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/chr${CHR}unmapped.bed

awk 'BEGIN {OFS="\t"} {print $3, $4, $5}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/genetic_map_chr${CHR}_combined_hg38_mapped.bed > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/tmpfile

## Sort position in increasing order
sort -k1,1n  /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/tmpfile > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/genetic_map_chr${CHR}_combined_hg38_mapped.txt

plink --bfile  /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/dyslypidemia_EUR_chr${CHR} --cm-map /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/genetic_map_chr@_combined_hg38_mapped.txt --make-bed --out $PWD/ldscore/dyslypidemia_EUR_chr${CHR}_ldscore

/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldsc/ldsc.py --bfile $PWD/ldscore/dyslypidemia_EUR_chr${CHR}_ldscore --l2 --ld-wind-cm 1 --yes-really --out $PWD/ldscore/dyslipidemia_EUR_chr${CHR}_LD
done


for CHR in {1..22}; do
zcat dyslipidemia_EUR_chr${CHR}_LD.l2.ldscore.gz |awk 'NR==1 {print $1, $2, $3, $4} NR>1 {sub("chr", "", $2); print $1 "\t" $2 "\t" $3 "\t" $4}' | gzip >  EUR/${CHR}.l2.ldscore.gz
cp dyslipidemia_EUR_chr${CHR}_LD.l2.M_5_50 EUR/${CHR}.l2.M_5_50
done




## Run MTAG
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/mtag.py --sumstats eur_chltot.shared.txt,eur_ldl.shared.txt,eur_hdl.shared.txt,eur_trigly.shared.txt,eur_nonhdl.shared.txt --snp_name snpid --chr_name chr --bpos_name bpos --z_name STAT --ld_ref_panel /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/EUR/ --out ./EUR_chltot_ldl_hdl_trigly_nonhdl_mtag --force --n_min 0.0 --stream_stdout

for trait in {1..5}; do
 sort -k12,12g EUR_chltot_ldl_hdl_trigly_nonhdl_mtag_trait_${trait}.txt > EUR_chltot_ldl_hdl_trigly_nonhdl_mtag_trait_${trait}_sorted.txt
done

#############
## African ##
#############
for CHR in {1..22}; do
plink --bfile  /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/dyslypidemia_AFR_chr${CHR} --cm-map /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/1000GP_Phase3/genetic_map_chr@_combined_hg38_mapped.txt --make-bed --out $PWD/ldscore/dyslypidemia_AFR_chr${CHR}_ldscore

/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldsc/ldsc.py --bfile $PWD/ldscore/dyslypidemia_AFR_chr${CHR}_ldscore --l2 --ld-wind-cm 1 --yes-really --out $PWD/ldscore/dyslipidemia_AFR_chr${CHR}_LD
done

for CHR in {1..22}; do
zcat dyslipidemia_AFR_chr${CHR}_LD.l2.ldscore.gz |awk 'NR==1 {print $1, $2, $3, $4} NR>1 {sub("chr", "", $2); print $1 "\t" $2 "\t" $3 "\t" $4}' | gzip >  AFR/${CHR}.l2.ldscore.gz
cp dyslipidemia_AFR_chr${CHR}_LD.l2.M_5_50 AFR/${CHR}.l2.M_5_50
done


## Run MTAG
/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/mtag/mtag.py --sumstats afr_chltot.shared.txt,afr_ldl.shared.txt,afr_hdl.shared.txt,afr_trigly.shared.txt,afr_nonhdl.shared.txt --snp_name snpid --chr_name chr --bpos_name bpos --z_name STAT --ld_ref_panel /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/ldscore/AFR/ --out ./AFR_chltot_ldl_hdl_trigly_nonhdl_mtag_2 --force --n_min 0.0 --stream_stdout


for trait in {1..5}; do
 sort -k12,12g AFR_chltot_ldl_hdl_trigly_nonhdl_mtag_2_trait_${trait}.txt > AFR_chltot_ldl_hdl_trigly_nonhdl_mtag_2_trait_${trait}_sorted.txt
done

https://github.com/JonJala/mtag/issues/117