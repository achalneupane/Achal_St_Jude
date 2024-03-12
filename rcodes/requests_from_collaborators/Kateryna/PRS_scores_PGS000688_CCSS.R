## Prepare files
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/")
library(data.table)
library(dplyr)
#######################
## PGS000688 _GRCh38 ##
#######################

PGS000688 <- read.table("PGS000688_hmPOS_GRCh38.txt", sep = "\t", header = T)
PGS000688.saved <- PGS000688[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight", "chr_position")]
colnames(PGS000688.saved) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "POS_GRCh37")

PGS000688 <- PGS000688[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight")]
PGS000688$TYPE <- "PGS000688"
PGS000688$TYPE2 <- "PGS000688"
colnames(PGS000688) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")


################################
BCC <- PGS000688

## Swap alleles if beta is negative
BCC[BCC$Effect_size < 0, c("REF", "Effect_allele")] <- BCC[BCC$Effect_size < 0, c("Effect_allele", "REF")]
BCC$Effect_size <- abs(BCC$Effect_size)

BCC$POS_GRCh38 <- as.numeric(BCC$POS_GRCh38)
BCC$CHROM <- paste0("chr", BCC$CHROM)
write.table(BCC, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/GRCh38_PGS000688.dat", col.names = F, row.names = F, quote = F)

BCC$KEY <- paste0(BCC$CHROM, ":", BCC$POS_GRCh38)

#############################
## Creat Key for Bim match ##
#############################
BCC$KEY <- paste0(BCC$CHROM, ":", BCC$POS_GRCh38)
BCC$REF_flipped <- chartr("acgtACGT", "tgcaTGCA", BCC$REF)
BCC$Effect_allele_flipped <- chartr("acgtACGT", "tgcaTGCA", BCC$Effect_allele)

BCC$KEY2 <- paste0(BCC$KEY,":", BCC$REF,":",BCC$Effect_allele, ";", BCC$KEY,":", BCC$Effect_allele, ":", BCC$REF, ";", BCC$KEY,":", BCC$REF_flipped,":",BCC$Effect_allele_flipped, ";", BCC$KEY,":", BCC$Effect_allele_flipped, ":", BCC$REF_flipped,";")


####################################
## extract variants from ccss_exp ## ## GRCh38
####################################

BIM.ALL <- {}
for (CHR in 1:22){
  print(paste0("Doing CHR: ", CHR))
  BIM <- fread(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr",CHR, "_ID_edited.bim"), header = F)
  BIM$KEY <- paste0("chr", BIM$V1, ":", BIM$V4)
  BIM <- BIM[BIM$KEY %in% BCC$KEY,]
  BIM.ALL <- rbind.data.frame(BIM.ALL,BIM)
}

# cc <- BCC[BCC$CHROM == "chr1",]
# cc <- cc [!duplicated(cc$KEY),]

BIM.final <- BIM.ALL[!duplicated(BIM.ALL$V2),]
## 0


# ## Use rowwise method below if too many variants. 
# BCC$BIM_CCSS_exp <- NA
# for (i in 1:nrow(BCC)){
#   snp_id <- BCC$KEY2[i]
#   snps <- unlist(strsplit(snp_id, ";"))
#   for (k in 1:length(snps)){
#     matching_row <- which(BIM.final$V2 %in% snps)
#   }
#   # If a match is found, assign the value to BCC$BIM
#   if (length(matching_row) > 0) {
#     BCC$BIM_CCSS_exp[i] <- BIM.final$V2[matching_row]
#   }
# }

cc <- BCC %>%
  rowwise() %>%
  mutate(BIM_CCSS_exp = paste(BIM.final$V2[BIM.final$V2 %in% unlist(strsplit(KEY2, ";"))], collapse = ";"))
cc$BIM_CCSS_exp[cc$BIM_CCSS_exp == ""] <- NA


# sum(cc$BIM_CCSS_exp == BCC$BIM_CCSS_exp, na.rm = T)
wanted.vars <- cc$BIM_CCSS_exp[!is.na(cc$BIM_CCSS_exp)]

## extract ccss_org variants
write.table(wanted.vars, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/extract_CCSS_exp_vars", col.names = F, row.names = F, quote = F)


####################################
## extract variants from ccss_org ## ## GRCh37
####################################
BCC.GRCh37 <- PGS000688.saved
BCC.GRCh37 <- BCC.GRCh37[c("CHROM", "POS_GRCh37", "REF", "Effect_allele", "Effect_size")]
BCC.GRCh37$TYPE <- "PGS000688"
BCC.GRCh37$TYPE2 <- "PGS000688"

BCC.GRCh37[BCC.GRCh37$Effect_size < 0, c("REF", "Effect_allele")] <- BCC.GRCh37[BCC.GRCh37$Effect_size < 0, c("Effect_allele", "REF")]
BCC.GRCh37$Effect_size <- abs(BCC.GRCh37$Effect_size)

BCC.GRCh37$POS_GRCh37 <- as.numeric(BCC.GRCh37$POS_GRCh37)
BCC.GRCh37$CHROM <- paste0("chr", BCC.GRCh37$CHROM)

write.table(BCC.GRCh37, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/GRCh37_PGS000688.dat", col.names = F, row.names = F, quote = F)

BCC.GRCh37$KEY <- paste0(BCC.GRCh37$CHROM, ":", BCC.GRCh37$POS_GRCh37)

BCC.GRCh37$REF_flipped <- chartr("acgtACGT", "tgcaTGCA", BCC.GRCh37$REF)
BCC.GRCh37$Effect_allele_flipped <- chartr("acgtACGT", "tgcaTGCA", BCC.GRCh37$Effect_allele)
BCC.GRCh37$KEY2 <- paste0(BCC.GRCh37$KEY,":", BCC.GRCh37$REF,":",BCC.GRCh37$Effect_allele, ";", BCC.GRCh37$KEY,":", BCC.GRCh37$Effect_allele, ":", BCC.GRCh37$REF, ";", BCC.GRCh37$KEY,":", BCC.GRCh37$REF_flipped,":",BCC.GRCh37$Effect_allele_flipped, ";", BCC.GRCh37$KEY,":", BCC.GRCh37$Effect_allele_flipped, ":", BCC.GRCh37$REF_flipped,";")


## First rename all BIM files with CHR:POS:REF:ALT in shell script in CCSS_org bim
BIM.ALL <- {}
for (CHR in 1:22){
  print(paste0("Doing CHR: ", CHR))
  BIM <- fread(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/CCSS_org_GRCh37_chr",CHR, ".bim"), header = F)
  BIM$KEY <- paste0(BIM$V1, ":", BIM$V4)
  BIM <- BIM[BIM$KEY %in% BCC.GRCh37$KEY,]
  BIM.ALL <- rbind.data.frame(BIM.ALL,BIM)
}

# cc <- BCC[BCC$CHROM == "chr1",]
# cc <- cc [!duplicated(cc$KEY),]

BIM.final <- BIM.ALL[!duplicated(BIM.ALL$V2),]
## 0

cc <- BCC.GRCh37 %>%
  rowwise() %>%
  mutate(BIM_CCSS_org = paste(BIM.final$V2[BIM.final$V2 %in% unlist(strsplit(KEY2, ";"))], collapse = ";"))
cc$BIM_CCSS_org[cc$BIM_CCSS_org == ""] <- NA


wanted.vars <- cc$BIM_CCSS_org[!is.na(cc$BIM_CCSS_org)]

## extract ccss_org variants
write.table(wanted.vars, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/extract_CCSS_org_vars", col.names = F, row.names = F, quote = F)



