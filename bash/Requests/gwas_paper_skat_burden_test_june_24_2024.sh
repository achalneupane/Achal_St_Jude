cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
awk -F'\t' 'NR==1 || ($210 ~ /Pathogenic|Likely_pathogenic/) && $210 !~ /Conflicting|Benign/' \
chr*.preQC_biallelic_renamed_ID_edited_gnomAD_clinvar_12_10_2023_FIELDS-simple4.txt > all_new_clinvar_P_LP.txt

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024
# ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/all_new_clinvar_P_LP.txt .
# cut -f3 all_new_clinvar_P_LP.txt| cut -d';' -f1 > extract_PLP.txt
# extract variants for EUR and AFR
module load plink/1.90b
for i in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${i}.PASS.decomposed --extract extract_PLP.txt --keep EURsamples.list --max-maf 0.05 --keep-allele-order --make-bed --out tmp_NFE_filtered_variants_EUR_chr${i}.dat
done

## EUR
ls -1 tmp_NFE_filtered_variants_EUR_chr*bed| sort -V | sed 's/\.bed//' > merge_list.txt
plink --merge-list merge_list.txt --keep-allele-order --make-bed --out tmp_EUR_merged_data_plink
module load plink2
plink2 --bfile tmp_EUR_merged_data_plink --freq counts --out tmp_EUR_merged_data_plink_MAC_EUR


## Extract EUR with maf 0.05, mac >=2 and with at least 2 variants
plink --bfile tmp_EUR_merged_data_plink --extract vars.clinvar.missense.loftee.eur.maf0.05.mac3.txt --keep-allele-order --make-bed --out EUR_final

# Now run gene based analysis. First prepare datasets
mkdir EUR
cd EUR

ln -s ../EUR_final.* .
ln -s ../extract-variants-final-EUR.GENE-SNP .
# 1- CREATE DIRECTORY FOR EACH CHR
for chr in {1..22}; do  mkdir chr${chr}; done


# 2 - MAKE RAW FILES FOR EACH CHR AND MOVE THEM INTO EACH DIRECTORY
CHR_TO_PROCESS="$(seq 1 22)"
THREADS=4
module load parallel
# BFILE="EOAD-2828-logistic-BETA"
BFILE="EUR_final"
parallel -j ${THREADS} plink --bfile ${BFILE} --keep-allele-order --allow-no-sex --recode A --output-chr MT --chr {} --out chr{}/${BFILE}-chr{} ::: ${CHR_TO_PROCESS}

# 3 - Extract annot info per chr and create gene-subsets
for file in $(ls *.GENE-SNP); do
INFILE=${file}
for chr in chr{1..22}; do
awk -v chr=${chr} '{if($1==chr) print $2,$3}' ${INFILE} > ${chr}/${chr}-${INFILE}
done
done



# 4 - Generate gen-sets using script get_gene.perl
set="EUR"    
BASE="extract-variants-final"
HOMEDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/EUR"
for chr in {1..22}; do \
cd chr${chr}; mkdir geneset-${set}; perl /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/rare_variants/get_gene.perl chr${chr}-${BASE}-${set}.GENE-SNP; mv *.gene geneset-${set}; cd ${HOMEDIR}; done

# Now Run SKAT-loop-OPTIMAL.r

#########
## AFR ##
#########
module load plink/1.90b
for i in {1..22}; do
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${i}.PASS.decomposed --extract extract_PLP.txt --keep AFRsamples.list --max-maf 0.05 --keep-allele-order --make-bed --out tmp_AFR_filtered_variants_AFR_chr${i}.dat
done
ls -1 tmp_AFR_filtered_variants_AFR_chr*bed| sort -V | sed 's/\.bed//' > AFRmerge_list.txt
plink --merge-list AFRmerge_list.txt --keep-allele-order --make-bed --out tmp_AFR_merged_data_plink
module load plink2
plink2 --bfile tmp_AFR_merged_data_plink --freq counts --out tmp_AFR_merged_data_plink_MAC_AFR


## Extract AFR with maf 0.05, mac >=2 and with at least 2 variants
plink --bfile tmp_AFR_merged_data_plink --extract vars.clinvar.missense.loftee.afr.maf0.05.mac3.txt --keep-allele-order --make-bed --out AFR_final

# Now run gene based analysis. First prepare datasets
mkdir AFR
cd AFR

ln -s ../AFR_final.* .
ln -s ../extract-variants-final-AFR.GENE-SNP .
# 1- CREATE DIRECTORY FOR EACH CHR
for chr in {1..22}; do  mkdir chr${chr}; done


# 2 - MAKE RAW FILES FOR EACH CHR AND MOVE THEM INTO EACH DIRECTORY
CHR_TO_PROCESS="$(seq 1 22)"
THREADS=4
module load parallel
# BFILE="EOAD-2828-logistic-BETA"
BFILE="AFR_final"
parallel -j ${THREADS} plink --bfile ${BFILE} --keep-allele-order --allow-no-sex --recode A --output-chr MT --chr {} --out chr{}/${BFILE}-chr{} ::: ${CHR_TO_PROCESS}

# 3 - Extract annot info per chr and create gene-subsets
for file in $(ls *.GENE-SNP); do
INFILE=${file}
for chr in chr{1..22}; do
awk -v chr=${chr} '{if($1==chr) print $2,$3}' ${INFILE} > ${chr}/${chr}-${INFILE}
done
done



# 4 - Generate gen-sets using script get_gene.perl
set="AFR"    
BASE="extract-variants-final"
HOMEDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/AFR"
for chr in {1..22}; do \
cd chr${chr}; mkdir geneset-${set}; perl /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/rare_variants/get_gene.perl chr${chr}-${BASE}-${set}.GENE-SNP; mv *.gene geneset-${set}; cd ${HOMEDIR}; done

# Now Run SKAT-loop-OPTIMAL.r




## Next, perform association analysis with echo:

plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed --extract chr16:25611595:C:T.txt --keep EURsamples.list --keep-allele-order --make-bed --out chr16:25611595:C:T_EUR.dat
plink --bfile chr16:25611595:C:T_EUR.dat --keep-allele-order --allow-no-sex --recode A --out chr16_25611595_C_T_EUR.dat_recodeA

plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed --extract chr16:25611595:C:T.txt --keep AFRsamples.list --keep-allele-order --make-bed --out chr16:25611595:C:T_AFR.dat
plink --bfile chr16:25611595:C:T_AFR.dat --keep-allele-order --allow-no-sex --recode A --out chr16_25611595_C_T_AFR.dat_recodeA