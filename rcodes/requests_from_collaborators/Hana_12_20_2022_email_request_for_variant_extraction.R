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

## Add Frequency
FREQ <- read.table("LMNA_ALL_FREQ_result.frq", header = T)
sum(FREQ$SNP %in% dd$SNP_ID)
dd$NORTHWESTERN_WGS_MAF <- FREQ$MAF[match(dd$SNP_ID, FREQ$SNP)]

dd$VQSR_PASS <- ifelse(dd$SNP_ID %in% LMNA_VQSR.pASS$V2, "YES", "NO")
dd <- dd[c("SNP_ID", "CHROM",  "POS",    "REF",    "ALT",    "VQSR_PASS", "JW1-12", "JW10-1", "JW11-6", "JW12-17", "JW13-2", "JW2-2", "JW3-4", "JW4-5", "JW5-5", "JW6", "JW7-7", "JW8-6", "JW9-12", "NORTHWESTERN_WGS_MAF")]


## Add annotation
LMNA.anno <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/HANA_Northwestern_request/LMNA_ANNOVAR.txt", sep = "\t", header = T)
LMNA.anno <- cbind.data.frame(ID = paste0(LMNA.anno$Otherinfo4,":", LMNA.anno$Otherinfo5,":",
       LMNA.anno$Otherinfo7,":", LMNA.anno$Otherinfo8), LMNA.anno)


sum(dd$SNP_ID %in% LMNA.anno$ID)
# 452
dd$SNP_ID[!dd$SNP_ID %in% LMNA.anno$ID]
# [1] "chr1:156084329:G:*" "chr1:156084335:T:*" "chr1:156084338:T:*"
# "chr1:156106643:A:*" "chr1:156109819:C:*" "chr1:156109820:A:*"
# "chr1:156109837:C:*" "chr1:156109839:C:*"

# dd$SNP_UPDATED <- dd$SNP_ID
# dd$SNP_UPDATED [grepl("chr1:156084329:G:*", dd$SNP_UPDATED)] <- "chr1:156084328:CG:C"
# dd$SNP_UPDATED [grepl("chr1:156084335:T:*", dd$SNP_UPDATED)] <- "chr1:156084334:GT:G"
# dd$SNP_UPDATED [grepl("chr1:156084338:T:*", dd$SNP_UPDATED)] <- "chr1:156084337:GT:G"
# # dd$SNP_UPDATED [grepl("chr1:156106643:A:*", dd$SNP_UPDATED)] <- "chr1:156106643:A"
# # dd$SNP_UPDATED [grepl("chr1:156109819:C:*", dd$SNP_UPDATED)] <- "chr1:156109819:C:"
# # dd$SNP_UPDATED [grepl("chr1:156109820:A:*", dd$SNP_UPDATED)] <- "chr1:156109820:A"
# # dd$SNP_UPDATED [grepl("chr1:156109837:C:*", dd$SNP_UPDATED)] <- "chr1:156109837:C:*"
# # dd$SNP_UPDATED [grepl("chr1:156109839:C:*", dd$SNP_UPDATED)] <- "chr1:156109839:C"

dd <- cbind.data.frame(dd, LMNA.anno[match(dd$SNP_ID, LMNA.anno$ID),])





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

## Add Frequency
FREQ <- read.table("EMD_ALL_FREQ_result.frq", header = T)
sum(FREQ$SNP %in% dd$SNP_ID)
dd$NORTHWESTERN_WGS_MAF <- FREQ$MAF[match(dd$SNP_ID, FREQ$SNP)]

dd$VQSR_PASS <- ifelse(dd$SNP_ID %in% EMD_VQSR.pASS$V2, "YES", "NO")
dd <- dd[c("SNP_ID", "CHROM",  "POS",    "REF",    "ALT",    "VQSR_PASS", "JW1-12", "JW10-1", "JW11-6", "JW12-17", "JW13-2", "JW2-2", "JW3-4", "JW4-5", "JW5-5", "JW6", "JW7-7", "JW8-6", "JW9-12", "NORTHWESTERN_WGS_MAF")]


## Add annotation
EMD.anno <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/HANA_Northwestern_request/EMD_ANNOVAR.txt", sep = "\t", header = T)
EMD.anno <- cbind.data.frame(ID = paste0(EMD.anno$Otherinfo4,":", EMD.anno$Otherinfo5,":",
                                          EMD.anno$Otherinfo7,":", EMD.anno$Otherinfo8), EMD.anno)


sum(dd$SNP_ID %in% EMD.anno$ID)
# 14
dd$SNP_ID[!dd$SNP_ID %in% EMD.anno$ID]

dd$CHROM - NA
dd$CHROM <- c("X")

dd <- cbind.data.frame(dd,EMD.anno[match(dd$SNP_ID, EMD.anno$ID),])

dd <- dd[c("SNP_ID","CHROM","POS","REF","ALT","VQSR_PASS","JW1-12","JW10-1",
     "JW11-6","JW12-17","JW13-2","JW2-2","JW3-4","JW4-5","JW5-5","JW6",
     "JW7-7","JW8-6","JW9-12","NORTHWESTERN_WGS_MAF","Func.refGene",
     "Gene.refGene","ExonicFunc.refGene","AAChange.refGene","X1000g2015aug_all",
     "gnomAD_genome_ALL","gnomAD_genome_AFR","gnomAD_genome_AMR","gnomAD_genome_ASJ",
     "gnomAD_genome_EAS","gnomAD_genome_FIN","gnomAD_genome_NFE","gnomAD_genome_OTH",
     "cosmic70","nci60","CLNSIG","REVEL")]

write.table(dd, "Hana_EMD_genotype.txt", sep = "\t", col.names = T, quote = F, row.names = F)


## on 03/27/2023, Hana requested to double check
LMNA.Pos <- read.table("LMNA_Pos_key.txt", sep = "\t", header = T)
LMNA.Pos$KEY <- paste(LMNA.Pos$CHROM, LMNA.Pos$POS, LMNA.Pos$REF, LMNA.Pos$ALT, sep = ":")
jw6 <- read.table("chr1.genotype_JW6.txt")
jw6$KEY <- gsub("chr","", paste(jw6$V1, jw6$V2, jw6$V3, jw6$V4, sep = ":"))

JW6.LMNA <- cbind.data.frame(LMNA.Pos,jw6[match(LMNA.Pos$KEY, jw6$KEY),])


write.table(JW6.LMNA, "JW6.LMNA.geno.txt", sep = "\t", quote = F, row.names = F, col.names = T)


jw10_1 <- read.table("chr1.genotype_JW10-1.txt")
jw10_1$KEY <- gsub("chr","", paste(jw10_1$V1, jw10_1$V2, jw10_1$V3, jw10_1$V4, sep = ":"))
jw10_1.LMNA <- cbind.data.frame(LMNA.Pos,jw10_1[match(LMNA.Pos$KEY, jw10_1$KEY),])


write.table(jw10_1.LMNA, "jw10_1.LMNA.geno.txt", sep = "\t", quote = F, row.names = F, col.names = T)
