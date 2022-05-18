## Extract variants from genotype data give the list of RSIDs
# Yadav's email dated: 5/17/2022
# Email subject: SNP data
# Details: Can you please extract the genotypes of the SNPs among SJLIFE survivors in the Excel file? The Patient list has ids as MRNs which you can use to get their corresponding sjlids and then get genotypes. Once you get their genotypes based on sjlids, please replace them back to MRNs. If you could provide actual genotypes (eg, AT, TT, AA), that would be great.

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/tasks/task_5_17_2022
ln -s ../../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr*.PASS.decomposed_geno.0.1_hwe.1e-10.* .

...{r}
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/tasks/task_5_17_2022")
## read annotated SJLIFE annotated VCF 
library(data.table)
## For Yadav
# list of variants; first create patterns to be searched in dbSNP
yadav.vars <- c("17179101", "12917", "3212227", "3091339", "776746", "2071487", "1138272", "1695", "2266637", "1728", "429358", "6265", "1611115", "2519152", "4680", "6277", "1800496", "1800443", "1799835", "1137070", "6296", "3746544", "10513112", "341", "70991108", "10760502", "10106", "11545078", "1805087", "1801394", "2236225", "1950902", "1801133", "1801131", "1051266", "1979277", "34489327", "1333049", "6922269", "2943634", "17465637", "501120", "17228212", "2815063", "6689879", "3557", "3816", "2677", "21387", "3055", "963", "7305", "11984041", "2023938", "4959130", "17580", "12124533", "13143308", "12932445")
# yadav.vars <- paste0("RS=", yadav.vars,sep = ";", collapse = "|")
yadav.vars <- paste0("RS=", yadav.vars,sep = ";", collapse = "|")
paste(paste("rs", yadav.vars,sep = ""), collapse = ",")
## Now look for these variants in UCSC browser by replacing , with \n: https://genome.ucsc.edu/cgi-bin/hgTables
# rs17179101,rs12917,rs3212227,rs3091339,rs776746,rs2071487,rs1138272,rs1695,rs2266637,rs1728,rs429358,rs6265,rs1611115,rs2519152,rs4680,rs6277,rs1800496,rs1800443,rs1799835,rs1137070,rs6296,rs3746544,rs10513112,rs341,rs70991108,rs10760502,rs10106,rs11545078,rs1805087,rs1801394,rs2236225,rs1950902,rs1801133,rs1801131,rs1051266,rs1979277,rs34489327,rs1333049,rs6922269,rs2943634,rs17465637,rs501120,rs17228212,rs2815063,rs6689879,rs3557,rs3816,rs2677,rs21387,rs3055,rs963,rs7305,rs11984041,rs2023938,rs4959130,rs17580,rs12124533,rs13143308,rs12932445
...

# The output file is then saved as yadav_rs_id_to_var_position.txt

## Alternatively, grep in dbSNP
# zgrep -a -m 1 '^#CHROM' /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz| cut -d$'\t' -f1-5 > yadav_rs_id_to_var_position.txt
# zgrep -a -E "RS=17179101;|RS=12917;|RS=3212227;|RS=3091339;|RS=776746;|RS=2071487;|RS=1138272;|RS=1695;|RS=2266637;|RS=1728;|RS=429358;|RS=6265;|RS=1611115;|RS=2519152;|RS=4680;|RS=6277;|RS=1800496;|RS=1800443;|RS=1799835;|RS=1137070;|RS=6296;|RS=3746544;|RS=10513112;|RS=341;|RS=70991108;|RS=10760502;|RS=10106;|RS=11545078;|RS=1805087;|RS=1801394;|RS=2236225;|RS=1950902;|RS=1801133;|RS=1801131;|RS=1051266;|RS=1979277;|RS=34489327;|RS=1333049;|RS=6922269;|RS=2943634;|RS=17465637;|RS=501120;|RS=17228212;|RS=2815063;|RS=6689879;|RS=3557;|RS=3816;|RS=2677;|RS=21387;|RS=3055;|RS=963;|RS=7305;|RS=11984041;|RS=2023938;|RS=4959130;|RS=17580;|RS=12124533;|RS=13143308;|RS=12932445;"  /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/snpEff/data/dbSNP/dbSNP_155.GCF_000001405.39.gz | cut -d$'\t' -f1-5 >> yadav_rs_id_to_var_position.txt

```{r}
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/tasks/task_5_17_2022")
df <- read.table("yadav_rs_id_to_var_position.txt", sep = "\t", header = T, comment.char = "")
dim(df)
head(df)

library(stringr)
df$X.chrom <- str_split(df$X.chrom, "_", simplify=T)[,1]
## Make SNPs biallelic
df$alts <- gsub("\\,$", "", df$alts)
library(tidyr)
df <- separate_rows(df,alts,sep=",")
# df$SNP <- gsub("chr", "", paste(df$X.chrom,df$chromEnd,df$ref, df$alts, sep = ":"))
# df$SNP_extract <- gsub("chr", "", paste(df$X.chrom,df$chromEnd, sep = ":"))
df$SNP <- paste(df$X.chrom,df$chromEnd,df$ref, df$alts, sep = ":")
df$SNP_extract <- paste(df$X.chrom,df$chromEnd, sep = ":")

## Extract these variants with PLINK (there could be allele flips, so extracting by chr and pos)
write.table(unique(df$SNP_extract), "/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/tasks/task_5_17_2022/SNPs_to_extract_05-18-2022-task.txt", row.names = F, col.names = F, quote = FALSE)
```

# Check which ones are multi allelic in our data
for line in $(cat SNPs_to_extract_05-18-2022-task.txt); do
echo "---- DOING ${line} ---"
grep -w ${line} MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr*.PASS.decomposed_geno.0.1_hwe.1e-10.bim
done

## This is multi allelic:
# ---- DOING chr2:226203364 ---
# MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr2.PASS.decomposed_geno.0.1_hwe.1e-10.bim:2     chr2:226203364:A:C      0       226203364       C       A
# MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr2.PASS.decomposed_geno.0.1_hwe.1e-10.bim:2     chr2:226203364:A:G      0       226203364       G       A



# Extract variants from plink
module load plink/1.90b
for i in {1..22}; do
BFILE="MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr${i}.PASS.decomposed_geno.0.1_hwe.1e-10"	
plink --bfile ${BFILE} --extract SNPs_to_extract_05-18-2022-task.txt --keep-allele-order --make-bed --out CHR_${i}
done