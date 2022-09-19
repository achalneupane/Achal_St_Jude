# Main folder: /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")
# load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")




splice.region <- LoF.unique[grepl("splice_region_variant", LoF.unique$`ANN[*].EFFECT`),]
dim(splice.region)
# 242541     25
non.splice.region <- LoF.unique[!grepl("splice_region_variant", LoF.unique$`ANN[*].EFFECT`),]
dim(non.splice.region)
# 73845    25

NO.splice.region <- splice.region[grepl("frameshift|stop|donor|acceptor", splice.region$`ANN[*].EFFECT`),]
dim(NO.splice.region)
# 4421   25

LOF.keep <- rbind.data.frame(non.splice.region, NO.splice.region)
dim(LOF.keep)
# 78266    25

# remove any missense that is not coming from clinvar
Predicted.vars.inVCF.Unique <-  rbind.data.frame(CLINVAR.unique, LOF.keep)

FINAL <- {}

CHROM <- 1:22
for (i in 1:length(CHROM)){
print(paste0("DOING chr", CHROM[i]))
## Extract gnomAD MAF
annovar.chr <- read.delim(paste0("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/gnomAD_maf_lt_0.01_chr",CHROM[i],".txt"), sep = "\t", header = T)
# annovar.chr <- read.delim(paste0("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01_chr",CHROM[i],".txt"), sep = "\t", stringsAsFactors = F, header = T)

annovar.chr$KEY <- paste(annovar.chr$Otherinfo4, annovar.chr$Otherinfo5, annovar.chr$Otherinfo7, annovar.chr$Otherinfo8, sep = ":")
# annovar.chr$KEY.chr.pos <- paste(annovar.chr$Otherinfo4, annovar.chr$Otherinfo5, sep = ":")
# Predicted.vars.inVCF.Unique$KEY.chr.pos <- paste(Predicted.vars.inVCF.Unique$CHROM, Predicted.vars.inVCF.Unique$POS, sep = ":")

sum(Predicted.vars.inVCF.Unique$KEY %in% annovar.chr$KEY)
# sum(Predicted.vars.inVCF.Unique$KEY.chr.pos %in% annovar.chr$KEY.chr.pos)

# gnomAD ALL
Predicted.vars.inVCF.Unique$gnomAD_genome_ALL <- annovar.chr$gnomAD_genome_ALL [match(Predicted.vars.inVCF.Unique$KEY, annovar.chr$KEY)]
Predicted.vars.inVCF.Unique$gnomAD_genome_NFE <- annovar.chr$gnomAD_genome_NFE [match(Predicted.vars.inVCF.Unique$KEY, annovar.chr$KEY)]
Predicted.vars.inVCF.Unique$gnomAD_genome_AFR <- annovar.chr$gnomAD_genome_AFR [match(Predicted.vars.inVCF.Unique$KEY, annovar.chr$KEY)]

wanted.vars.all <- Predicted.vars.inVCF.Unique[!is.na(Predicted.vars.inVCF.Unique$gnomAD_genome_ALL),]
wanted.vars.all <- wanted.vars.all[wanted.vars.all$gnomAD_genome_ALL < 0.01,]

FINAL <- rbind.data.frame(FINAL, wanted.vars.all)
print(paste0("NROWS of FINAL ", nrow(FINAL)))
}
rm (annovar.chr)

# Keep the unique variants only
FINAL <- FINAL[!duplicated(FINAL$KEY),]


# Get the BIM ID
FINAL$BIM_ID <- NA

library(data.table)

CHROM <- 1:22
for (i in 1:length(CHROM)){
BIM <- fread(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr", CHROM[i], ".PASS.decomposed_geno.0.1_hwe.1e-10.bim"))
BIM$KEY.chr.pos <- paste0("chr",BIM$V1,":",BIM$V4)
BIM$KEY1 <- paste0("chr",BIM$V1,":",BIM$V4,":",BIM$V6,":", BIM$V5)
BIM$KEY2 <- paste0("chr",BIM$V1,":",BIM$V4,":",BIM$V5,":", BIM$V6)
sum(FINAL$KEY %in% BIM$KEY1) #2767
sum(FINAL$KEY %in% BIM$KEY2) # 22
sum(FINAL$KEY.chr.pos %in% BIM$KEY.chr.pos) # 2791

# FINAL$exact_match <- ifelse(FINAL$KEY %in% BIM$KEY1, "YES", "NO")
# FINAL$swapped_match <- ifelse(FINAL$KEY %in% BIM$KEY2, "YES", "NO")

# FINAL$exact_match_var <- BIM$KEY1[match(FINAL$KEY, BIM$KEY1)]
# FINAL$swapped_match_var <- BIM$KEY2[match(FINAL$KEY, BIM$KEY2)]

FINAL$BIM_ID[is.na(FINAL$BIM_ID)] <- BIM$V2 [match(FINAL$KEY[is.na(FINAL$BIM_ID)], BIM$KEY1)]
FINAL$BIM_ID[is.na(FINAL$BIM_ID)] <- BIM$V2 [match(FINAL$KEY [is.na(FINAL$BIM_ID)], BIM$KEY2)]
}

dim(FINAL)
# 38293    29
gene.table <- as.data.frame(table(FINAL$`ANN[*].GENE`))
## Keep only those with at least two variants
gene.table <- gene.table[gene.table$Freq >= 2,]

sum(FINAL$`ANN[*].GENE` %in% gene.table$Var1)
# 30709

save.original <- FINAL

saveRDS(FINAL, "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract_variants_gnomAD_ALL_NFE_lt_0.01.rds")

# FINAL <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract_variants_gnomAD_ALL_NFE_lt_0.01.rds")

FINAL <- FINAL[!is.na(FINAL$BIM_ID),]


## Extract NFE and AFR and get MAC for each variant from NFE and AFR datasets
FINAL.NFE <- FINAL[FINAL$gnomAD_genome_NFE < 0.01,]
FINAL.AFR <- FINAL[FINAL$gnomAD_genome_AFR < 0.01,]


write.table(as.data.frame(FINAL.NFE$BIM_ID), "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract_variants_gnomAD_ALL_NFE_lt_0.01_vars.list", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(as.data.frame(FINAL.AFR$BIM_ID), "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract_variants_gnomAD_ALL_AFR_lt_0.01_vars.list", sep = "\t", row.names = F, quote = F, col.names = F)



## step 2. 
## Get MAC

# test_MAC <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/test_dat.raw", stringsAsFactors = F, header = T)
# colSums(test_MAC[7:16], na.rm = T)

EUR_MAC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/MAC_EUR.acount", header = T, stringsAsFactors = F, skip = 'CHROM')
dim(EUR_MAC)
# 21272     6

AFR_MAC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/MAC_AFR.acount", header = T, stringsAsFactors = F, skip = 'CHROM')
dim(AFR_MAC)
# 26910     6

## Selecting variants that have MAC of at least 3
EUR_MAC <- EUR_MAC[EUR_MAC$ALT_CTS >= 3,]
dim(EUR_MAC)
# 4784    6
AFR_MAC <- AFR_MAC[AFR_MAC$ALT_CTS >= 3,]
dim(AFR_MAC)
# 2471    6


## Keeping only those variants with MAC of 3 or more
FINAL.NFE <- FINAL.NFE[FINAL.NFE$BIM_ID %in% EUR_MAC$ID,]
FINAL.AFR <- FINAL.AFR[FINAL.AFR$BIM_ID %in% AFR_MAC$ID,]

#############################################################
## Keeping only those with 2 or more variants in each gene ##
#############################################################
library(dplyr)

## NFE
FINAL.NFE <- FINAL.NFE[!duplicated(FINAL.NFE$BIM_ID),]
gene.counts.NFE <-   FINAL.NFE %>% 
group_by(`ANN[*].GENE`) %>% 
  count(`ANN[*].GENE`) 

FINAL.NFE$gene.counts <- gene.counts.NFE$n[match(FINAL.NFE$`ANN[*].GENE`, gene.counts.NFE$`ANN[*].GENE`)]

FINAL.NFE <- FINAL.NFE[FINAL.NFE$gene.counts >= 2,]
dim(FINAL.NFE)
# 1593   35

## AFR
FINAL.AFR <- FINAL.AFR[!duplicated(FINAL.AFR$BIM_ID),]
gene.counts.AFR <-   FINAL.AFR %>% 
  group_by(`ANN[*].GENE`) %>% 
  count(`ANN[*].GENE`) 

FINAL.AFR$gene.counts <- gene.counts.AFR$n[match(FINAL.AFR$`ANN[*].GENE`, gene.counts.AFR$`ANN[*].GENE`)]

FINAL.AFR <- FINAL.AFR[FINAL.AFR$gene.counts >= 2,]
dim(FINAL.AFR)
# 481  35

## This is the final file for analysis
## Gene list
EUR <- cbind.data.frame(CHROM=FINAL.NFE$CHROM, GENE=FINAL.NFE$`ANN[*].GENE`, SNPID=FINAL.NFE$BIM_ID)


AFR <- cbind.data.frame(CHROM=FINAL.AFR$CHROM, GENE=FINAL.AFR$`ANN[*].GENE`, SNPID=FINAL.AFR$BIM_ID)


write.table(EUR, "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract-variants-final-EUR.GENE-SNP", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(AFR, "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract-variants-final-AFR.GENE-SNP", sep = "\t", row.names = F, quote = F, col.names = F)

## FINAL VARS 
write.table(as.data.frame(EUR$SNPID), "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract_variants_gnomAD_ALL_NFE_lt_0.01_MAC_gt_3_at_least_2_vars_per_gene.list", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(as.data.frame(AFR$SNPID), "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract_variants_gnomAD_ALL_AFR_lt_0.01_MAC_gt_3_at_least_2_vars_per_gene.list", sep = "\t", row.names = F, quote = F, col.names = F)



# # Fisher test of AFR_diabetes_chrALL.dat-chr11.raw for chr11:71918924:G:A and chr11:71918995:TTTAGCAAGTTCTCTACAGGA:T
# df.raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/gene-based-analysis/AFR/chr11/AFR_diabetes_chrALL.dat-chr11.raw", header = T)
# df.raw <- grep("71918924", colnames(df.raw))

