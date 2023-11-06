cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/1kGP
module load plink/1.90b
plink --bfile 1000genomes_merged --keep AFR1KG --make-bed --out 1000genomes_AFR

# Step 1: Calculate LD with the reference SNP (rs6689879)
plink --bfile 1000genomes_AFR --r2 --ld-snp rs6689879 --ld-window-kb 100 --ld-window 99999 --ld-window-r2 0.6 --out ld_calculations

# Step 2: Extract SNPs in LD and meet MAF criteria
plink --bfile 1000genomes_AFR --extract ld_calculations.ld --maf 0.01 --make-bed --out extracted_snps

# Extract annotation for these variants
head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt > extracted_annotation_yadav.txt
egrep "$(awk '{print $2}' extracted_snps.bim | awk '{ORS = "|"; print}' | sed 's/|$//')" /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt >> extracted_annotation_yadav.txt


bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/1kGP//extracted_snps.bim")
df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/1kGP/extracted_annotation_yadav.txt", sep = "\t", header = T)
df <- df[c(1,2,3,4,5,7,8,9,11)]
df$rsid <- sapply(strsplit(df$ID, ";"), function(x) tail(unlist(x), 1))
df$ID <- sub(";rs.*$", "", df$ID)
df.unique <- df[!duplicated(df$rsid),]

sum(df.unique$rsid %in% bim$V2)
# 19

df.unique <- df.unique[df.unique$rsid %in% bim$V2,]
write.table(df.unique, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/1kGP/final_extracted_annotation_yadav.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# Thanks! Can you please expand the search to +/- 250 kb and also provide the following for each SNPs?

# LD r2 value with respect to rs6689879.
# Distance from rs6689879.
# Alleles
# MAF

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/1kGP
module load plink/1.90b
plink --bfile 1000genomes_merged --keep AFR1KG --make-bed --out 1000genomes_AFR
plink --bfile 1000genomes_AFR --freq --out 1000genomes_AFR_freq

# Step 1: Calculate LD with the reference SNP (rs6689879)
plink --bfile 1000genomes_AFR --r2 --ld-snp rs6689879 --ld-window-kb 250 --ld-window 99999 --ld-window-r2 0.6 --out ld_calculations_250kb

# Step 2: Extract SNPs in LD and meet MAF criteria
plink --bfile 1000genomes_AFR --extract ld_calculations_250kb.ld --maf 0.01 --make-bed --out extracted_snps_250kb

# Extract annotation for these variants
head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt > extracted_250kb_annotation_yadav.txt
egrep "$(awk '{print $2}' extracted_snps_250kb.bim | awk '{ORS = "|"; print}' | sed 's/|$//')" /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt >> extracted_250kb_annotation_yadav.txt


egrep "$(awk '{print $2}' extracted_snps_250kb.bim | awk '{ORS = "|"; print}' | sed 's/|$//')" 1000genomes_AFR_freq.frq > 1000genomes_AFR_freq_250kb.frq


bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/1kGP//extracted_snps_250kb.bim")
df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/1kGP/extracted_250kb_annotation_yadav.txt", sep = "\t", header = T)
df <- df[c(1,2,3,4,5,7,8,9,11)]
df$rsid <- sapply(strsplit(df$ID, ";"), function(x) tail(unlist(x), 1))
df$ID <- sub(";rs.*$", "", df$ID)
df.unique <- df[!duplicated(df$rsid),]

sum(df.unique$rsid %in% bim$V2)
# 19

df.unique <- df.unique[df.unique$rsid %in% bim$V2,]

## add freq
freq <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/1kGP/1000genomes_AFR_freq_250kb.frq")
df.unique$MAF <- freq$V5[match(df.unique$rsid , freq$V2)]

## add ld
ldscore <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/1kGP/ld_calculations_250kb.ld", header = T)
df.unique$LD <- ldscore$R2[match(df.unique$rsid , ldscore$SNP_B)]

df.unique$distance_BP <- 113559097 - df.unique$POS

write.table(df.unique, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/1kGP/final_extracted_250kb_with_annotation_yadav.txt", sep = "\t", col.names = T, row.names = F, quote = F)
