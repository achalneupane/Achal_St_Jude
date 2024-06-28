Step 1: Install PRSice
Install R: PRSice requires R to be installed on your system. You can download and install R from CRAN.

Install PRSice:

Download the PRSice software from the PRSice website.
Extract the downloaded zip file to a directory of your choice.
Step 2: Prepare Your Data
Base GWAS Summary Statistics: Obtain the summary statistics from a Genome-Wide Association Study (GWAS). The file should contain the following columns (column names may vary):

SNP (or rsID): The SNP identifier.
CHR: Chromosome number.
BP: Base pair position.
A1: Allele 1.
A2: Allele 2.
OR/BETA: Effect size (odds ratio or beta value).
P: p-value of the association.
Target Sample Genotype Data: Prepare your target sample genotype data in PLINK binary format (BED, BIM, and FAM files).

Phenotype File: A phenotype file containing the IDs and corresponding phenotypes. The file should be in the format required by PRSice.

Step 3: Running PRSice
Open Terminal/Command Prompt: Navigate to the directory where you extracted PRSice.

Command Syntax: Use the following command to run PRSice:




cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia
ln -s ../../PGS000699_hmPOS_GRCh38.txt.gz .
ln -s ../../PGS000686_hmPOS_GRCh38.txt.gz .
ln -s ../../PGS000688_hmPOS_GRCh38.txt.gz .
ln -s ../../PGS000888_hmPOS_GRCh38.txt.gz .
ln -s ../../PGS000677_hmPOS_GRCh38.txt.gz .

zcat PGS000699_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000699_hmPOS_GRCh38.txt
zcat PGS000686_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000686_hmPOS_GRCh38.txt
zcat PGS000688_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000688_hmPOS_GRCh38.txt
zcat PGS000888_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000888_hmPOS_GRCh38.txt
zcat PGS000677_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000677_hmPOS_GRCh38.txt

# Then run prsice_file_preparation.R


#!/bin/bash

# # Use awk to process the input file and create the new column
# awk '
# BEGIN {
#     FS=OFS="\t"
# }
# NR == 1 {
#     $NF=$NF"\tSNP\tP"
#     print
#     next
# }
# {
#     $NF=$NF"\t"$2":"$3":"$5":"$4"\t"
#     print
# }
# ' "PGS000677_hmPOS_GRCh38.txt" > "PGS000677_hmPOS_GRCh38_processed.txt"




# awk '
# BEGIN {
#     FS=OFS="\t"
# }
# NR == 1 {
#     $NF=$NF"\tSNP\tP"
#     print
#     next
# }
# {
#     $NF=$NF"\t"$8":"$9":"$4":"$3"\t"
#     print
# }
# ' "PGS000888_hmPOS_GRCh38.txt" > "PGS000888_hmPOS_GRCh38_processed.txt"

# split by chrmosomes
ls /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr*_ID_edited.bed|sort -V | sed 's/\.bed//' > target_list.txt

module load R/4.3.2-shlib

for file in PGS000699 PGS000686 PGS000688 PGS000888 PGS000677; do
echo "Doing:: ${file}"
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/prsice/PRSice_linux/PRSice.R \
--dir /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_exp/ \
--prsice /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/prsice/PRSice_linux/PRSice_linux \
--a1 effect_allele \
--a2 other_allele \
--chr CHR \
--bp BP \
--snp SNP \
--base ${file}_hmPOS_GRCh38_processed.txt \
--target-list target_list.txt \
--no-clump \
--no-regress \
--keep-ambig \
--stat BETA \
--pvalue P \
--fastscore \
--bar-levels 1 \
--type bed \
--out ${file}_PRS_output
done
# Rscript PRSice.R --dir . --prsice ./PRSice_linux --base TOY_BASE_GWAS.assoc --target TOY_TARGET_DATA --binary-target F --pheno TOY_TARGET_DATA.pheno --pheno-col Pheno --stat OR --type bed --out PRS_output

# https://www.biostars.org/p/9463113/
# https://groups.google.com/g/prsice/c/hzOugnwTmek
# https://www.linkedin.com/in/shing-wan-choi-13614aa5/

###################
## CCSS original ##
###################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_org

ln -s ../../PGS000699_hmPOS_GRCh38.txt.gz .
ln -s ../../PGS000686_hmPOS_GRCh38.txt.gz .
ln -s ../../PGS000688_hmPOS_GRCh38.txt.gz .
ln -s ../../PGS000888_hmPOS_GRCh38.txt.gz .
ln -s ../../PGS000677_hmPOS_GRCh38.txt.gz .

zcat PGS000699_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000699_hmPOS_GRCh38.txt
zcat PGS000686_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000686_hmPOS_GRCh38.txt
zcat PGS000688_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000688_hmPOS_GRCh38.txt
zcat PGS000888_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000888_hmPOS_GRCh38.txt
zcat PGS000677_hmPOS_GRCh38.txt.gz | grep -v "#" > PGS000677_hmPOS_GRCh38.txt


# Then run prsice_file_preparation.R


# split by chrmosomes
ls /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/CCSS_org_GRCh37_chr*.bed|sort -V | sed 's/\.bed//' > target_list.txt

module load R/4.3.2-shlib

for file in PGS000699 PGS000686 PGS000688 PGS000888 PGS000677; do
echo "Doing:: ${file}"
Rscript /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/prsice/PRSice_linux/PRSice.R \
--dir /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_org/ \
--prsice /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/prsice/PRSice_linux/PRSice_linux \
--a1 effect_allele \
--a2 other_allele \
--chr CHR \
--bp BP \
--snp SNP \
--base ${file}_hmPOS_GRCh38_processed.txt \
--target-list target_list.txt \
--no-clump \
--no-regress \
--keep-ambig \
--stat BETA \
--pvalue P \
--fastscore \
--bar-levels 1 \
--type bed \
--out ${file}_PRS_output
done