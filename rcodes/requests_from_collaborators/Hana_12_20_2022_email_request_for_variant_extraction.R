##################
## CCSS_exp_WGS ##
##################

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/HANA_Northwestern_request")

library(dplyr)

##########
## LMNA ##
##########
raw <- read.delim("LMNA_ALL_recodeA.raw", sep = " ", header = T, check.names = F)
dim(raw)
head(raw)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

# exclude chr9.84286010.A.G
raw <- raw[!grepl("FID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]
colnames(raw) # check this with the .bim REF and NON-REFERENCE alleles

rownames(raw) <- raw$IID
raw <- raw[!grepl("IID", colnames(raw))]

df.list <- list()
for (i in 1:ncol(raw)){
  # i=1
df.list.tmp <- raw[i]

REF = unlist(strsplit(colnames(df.list.tmp), "\\:"))[3]
ALT = unlist(strsplit(colnames(df.list.tmp), "\\:"))[4]

df.list.tmp$genotype <- as.character(df.list.tmp[,1])
df.list.tmp$genotype <- gsub("0", paste(REF,REF, sep = "/"), df.list.tmp$genotype)
df.list.tmp$genotype <- gsub("1", paste(REF,ALT, sep = "/"), df.list.tmp$genotype)
df.list.tmp$genotype <- gsub("2", paste(ALT,ALT, sep = "/"), df.list.tmp$genotype)
df.list.tmp$samples <- row.names(df.list.tmp)
colnames(df.list.tmp)
df.list.tmp[1] <- df.list.tmp$genotype
df.list.tmp <- df.list.tmp[-2]
df.list[[i]] <- df.list.tmp
}

cc <- Reduce(function(...) merge(..., by= "samples", all.x=T), df.list)
dd <- t(cc)
colnames(dd) <- dd[1,]
dd <- dd[-1,]
dd <- cbind.data.frame(SNP_ID = rownames(dd), dd)

LMNA.BIM <- read.table("LMNA_ALL.bim")
dd$CHROM <- LMNA.BIM$V1[match(dd$SNP_ID, LMNA.BIM$V2)]
dd$POS <- LMNA.BIM$V4[match(dd$SNP_ID, LMNA.BIM$V2)]
dd$REF <- LMNA.BIM$V6[match(dd$SNP_ID, LMNA.BIM$V2)]
dd$ALT <- LMNA.BIM$V5[match(dd$SNP_ID, LMNA.BIM$V2)]

# Flag VQSR PASS
LMNA_VQSR.pASS <- read.table("LMNA_ALL_VQSR_PASS.bim")

sum(dd$SNP_ID %in%   LMNA_VQSR.pASS$V2)
# 436

dd$VQSR_PASS <- ifelse(dd$SNP_ID %in% LMNA_VQSR.pASS$V2, "YES", "NO")
dd <- dd[c("CHROM",  "POS",    "REF",    "ALT",    "VQSR_PASS", "JW1-12", "JW10-1", "JW11-6", "JW12-17", "JW13-2", "JW2-2", "JW3-4", "JW4-5", "JW5-5", "JW6", "JW7-7", "JW8-6", "JW9-12")]

write.table(dd, "Hana_LMNA_genotype.txt", sep = "\t", col.names = T, quote = F, row.names = F)
#########
## EMD ##
#########

raw <- read.delim("EMD_ALL_recodeA.raw", sep = " ", header = T, check.names = F)
dim(raw)
head(raw)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

raw <- raw[!grepl("FID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]
colnames(raw) # check this with the .bim REF and NON-REFERENCE alleles

rownames(raw) <- raw$IID
raw <- raw[!grepl("IID", colnames(raw))]

df.list <- list()
for (i in 1:ncol(raw)){
  # i=1
  df.list.tmp <- raw[i]
  
  REF = unlist(strsplit(colnames(df.list.tmp), "\\:"))[3]
  ALT = unlist(strsplit(colnames(df.list.tmp), "\\:"))[4]
  
  df.list.tmp$genotype <- as.character(df.list.tmp[,1])
  df.list.tmp$genotype <- gsub("0", paste(REF,REF, sep = "/"), df.list.tmp$genotype)
  df.list.tmp$genotype <- gsub("1", paste(REF,ALT, sep = "/"), df.list.tmp$genotype)
  df.list.tmp$genotype <- gsub("2", paste(ALT,ALT, sep = "/"), df.list.tmp$genotype)
  df.list.tmp$samples <- row.names(df.list.tmp)
  colnames(df.list.tmp)
  df.list.tmp[1] <- df.list.tmp$genotype
  df.list.tmp <- df.list.tmp[-2]
  df.list[[i]] <- df.list.tmp
}

cc <- Reduce(function(...) merge(..., by= "samples", all.x=T), df.list)
dd <- t(cc)
colnames(dd) <- dd[1,]
dd <- dd[-1,]
dd <- cbind.data.frame(SNP_ID = rownames(dd), dd)

EMD.BIM <- read.table("EMD_ALL.bim")
dd$CHROM <- EMD.BIM$V1[match(dd$SNP_ID, EMD.BIM$V2)]
dd$POS <- EMD.BIM$V4[match(dd$SNP_ID, EMD.BIM$V2)]
dd$REF <- EMD.BIM$V6[match(dd$SNP_ID, EMD.BIM$V2)]
dd$ALT <- EMD.BIM$V5[match(dd$SNP_ID, EMD.BIM$V2)]

# Flag VQSR PASS
EMD_VQSR.pASS <- read.table("EMD_ALL_VQSR_PASS.bim")

sum(dd$SNP_ID %in%   EMD_VQSR.pASS$V2)
# 14

dd$VQSR_PASS <- ifelse(dd$SNP_ID %in% EMD_VQSR.pASS$V2, "YES", "NO")
dd <- dd[c("CHROM",  "POS",    "REF",    "ALT",    "VQSR_PASS", "JW1-12", "JW10-1", "JW11-6", "JW12-17", "JW13-2", "JW2-2", "JW3-4", "JW4-5", "JW5-5", "JW6", "JW7-7", "JW8-6", "JW9-12")]
dd$CHROM - NA
dd$CHROM <- c("X")

write.table(dd, "Hana_EMD_genotype.txt", sep = "\t", col.names = T, quote = F, row.names = F)
