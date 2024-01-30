## Prepare files
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/")
#########
## ST6 ##
#########
ST6 <- read.table("Choquet_ST6_multi_supplementary_6.txt", header = T, sep = "\t")

ST6.hg38 <- read.table("hgTables_ST6_multi.txt")
ST6$hg38BP <- ST6.hg38$V3[match(ST6$SNP, ST6.hg38$V4)]

ST6$hg38BP[ST6$SNP == "rs61824911"] <- "228843990"
ST6$hg38BP[ST6$SNP == "rs17401449"] <- "33480139"
ST6$hg38BP[ST6$SNP == "rs55804368"] <- "36630339"

ST6$beta <- log(ST6$OR)

ST6.saved <- ST6[,c("CHR", "hg38BP", "OA", "EA", "beta", "BP")]
colnames(ST6.saved) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "POS_GRCh37")
## Wanted columns
ST6 <- ST6[,c("CHR", "hg38BP", "OA", "EA", "beta")]
ST6$TYPE <- "ST6"
ST6$TYPE2 <- "ST6"
colnames(ST6) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")

###############
## PGS003416 ##
###############

PGS003416 <- read.table("PGS003416_hmPOS_GRCh38.txt", sep = "\t", header = T)
PGS003416.saved <- PGS003416[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight", "chr_position")]
colnames(PGS003416.saved) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "POS_GRCh37")

PGS003416 <- PGS003416[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight")]


PGS003416$TYPE <- "PGS003416"
PGS003416$TYPE2 <- "PGS003416"
colnames(PGS003416) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")

###############
## PGS000454 ##
###############

PGS000454 <- read.table("PGS000454_hmPOS_GRCh38.txt", sep = "\t", header = T)
PGS000454.saved <- PGS000454[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight", "chr_position")]
colnames(PGS000454.saved) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "POS_GRCh37")

PGS000454 <- PGS000454[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight")]
PGS000454$TYPE <- "PGS000454"
PGS000454$TYPE2 <- "PGS000454"
colnames(PGS000454) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")

###############
## PGS000356 ##
###############

PGS000356 <- read.table("PGS000356_hmPOS_GRCh38.txt", sep = "\t", header = T)
PGS000356.saved <- PGS000356[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight", "chr_position")]
colnames(PGS000356.saved) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "POS_GRCh37")

PGS000356 <- PGS000356[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight")]
PGS000356$TYPE <- "PGS000356"
PGS000356$TYPE2 <- "PGS000356"
colnames(PGS000356) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")


################################
BCC <- rbind.data.frame(ST6, PGS003416, PGS000454, PGS000356)

## Swap alleles if beta is negative
BCC[BCC$Effect_size < 0, c("REF", "Effect_allele")] <- BCC[BCC$Effect_size < 0, c("Effect_allele", "REF")]
BCC$Effect_size <- abs(BCC$Effect_size)

BCC$POS_GRCh38 <- as.numeric(BCC$POS_GRCh38)
BCC$CHROM <- paste0("chr", BCC$CHROM)
write.table(BCC, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/all_BCC_effect.dat", col.names = F, row.names = F, quote = F)

BCC$KEY <- paste0(BCC$CHROM, ":", BCC$POS_GRCh38)




## Swap alleles if beta is negative (GRCh37)
ST6.saved$TYPE <- "ST6"
ST6.saved$TYPE2 <- "ST6"

PGS003416.saved$TYPE <- "PGS003416"
PGS003416.saved$TYPE2 <- "PGS003416"

PGS000454.saved$TYPE <- "PGS000454"
PGS000454.saved$TYPE2 <- "PGS000454"

PGS000356.saved$TYPE <- "PGS000356"
PGS000356.saved$TYPE2 <- "PGS000356"

BCC.GRCh37 <- rbind.data.frame(ST6.saved, PGS003416.saved, PGS000454.saved, PGS000356.saved)

BCC.GRCh37[BCC.GRCh37$Effect_size < 0, c("REF", "Effect_allele")] <- BCC.GRCh37[BCC.GRCh37$Effect_size < 0, c("Effect_allele", "REF")]
BCC.GRCh37$Effect_size <- abs(BCC.GRCh37$Effect_size)

BCC.GRCh37$POS_GRCh37 <- as.numeric(BCC.GRCh37$POS_GRCh37)
BCC.GRCh37$CHROM <- paste0("chr", BCC.GRCh37$CHROM)

write.table(BCC.GRCh37, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/all_BCC_effect_ccss_org.dat", col.names = F, row.names = F, quote = F)

BCC.GRCh37$KEY <- paste0(BCC.GRCh37$CHROM, ":", BCC.GRCh37$POS_GRCh38)

#############################
## Creat Key for Bim match ##
#############################
BCC$KEY <- paste0(BCC$CHROM, ":", BCC$POS_GRCh38, ":")
BCC$REF_flipped <- chartr("acgtACGT", "tgcaTGCA", BCC$REF)
BCC$Effect_allele_flipped <- chartr("acgtACGT", "tgcaTGCA", BCC$Effect_allele)

BCC$KEY2 <- paste0(BCC$KEY,":", BCC$REF,":",BCC$Effect_allele, ";", BCC$KEY,":", BCC$Effect_allele, ":", BCC$REF, ";", BCC$KEY,":", BCC$REF_flipped,":",BCC$Effect_allele_flipped, ";", BCC$KEY,":", BCC$Effect_allele_flipped, ":", BCC$REF_flipped,";")


##################################
## extract variants from SJLIFE ##
##################################
library(data.table)
# BED <- cbind.data.frame(BCC$CHROM, BCC$POS_GRCh38-2, BCC$POS_GRCh38+2)
# write.table(BED, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/all_bed_BCC.bed", col.names = F, row.names = F, quote = F)

BIM.ALL <- {}
for (CHR in 1:22){
print(paste0("Doing CHR: ", CHR))
# BIM <- read.table(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr",CHR, ".PASS.decomposed.bim"), header = F)
BIM <- fread(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr",CHR, ".preQC_biallelic_renamed_ID_edited.vcf.gz.bim"), header = F)
BIM$KEY <- paste0("chr", BIM$V1, ":", BIM$V4)
BIM <- BIM[BIM$KEY %in% BCC$KEY,]
BIM.ALL <- rbind.data.frame(BIM.ALL,BIM)
}

# cc <- BCC[BCC$CHROM == "chr1",]
# cc <- cc [!duplicated(cc$KEY),]

BIM.final <- BIM.ALL[!duplicated(BIM.ALL$V2),]
## 0

BCC$BIM_SJLIFE <- NA
for (i in 1:nrow(BCC)){
  snp_id <- BCC$KEY2[i]
  snps <- unlist(strsplit(snp_id, ";"))
  for (k in 1:length(snps)){
  matching_row <- which(BIM.final$V2 %in% snps)
  }
  # If a match is found, assign the value to BCC$BIM
  if (length(matching_row) > 0) {
    BCC$BIM_SJLIFE[i] <- BIM.final$V2[matching_row]
  }
}


####################################
## extract variants from ccss_exp ##
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

BCC$BIM_CCSS_exp <- NA
for (i in 1:nrow(BCC)){
  snp_id <- BCC$KEY2[i]
  snps <- unlist(strsplit(snp_id, ";"))
  for (k in 1:length(snps)){
    matching_row <- which(BIM.final$V2 %in% snps)
  }
  # If a match is found, assign the value to BCC$BIM
  if (length(matching_row) > 0) {
    BCC$BIM_CCSS_exp[i] <- BIM.final$V2[matching_row]
  }
}

## extract ccss_org variants
write.table(BCC$BIM_SJLIFE[!is.na(BCC$BIM_SJLIFE)], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/extract_sjlife_vars", col.names = F, row.names = F, quote = F)

write.table(BCC$BIM_CCSS_exp[!is.na(BCC$BIM_CCSS_exp)], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/extract_CCSS_exp_vars", col.names = F, row.names = F, quote = F)


####################################
## extract variants from ccss_org ##
####################################
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

BCC.GRCh37$BIM_CCSS_org <- NA
for (i in 1:nrow(BCC.GRCh37)){
  snp_id <- BCC.GRCh37$KEY2[i]
  snps <- unlist(strsplit(snp_id, ";"))
  for (k in 1:length(snps)){
    matching_row <- which(BIM.final$V2 %in% snps)
  }
  # If a match is found, assign the value to BCC.GRCh37$BIM
  if (length(matching_row) > 0) {
    BCC.GRCh37$BIM_CCSS_org[i] <- BIM.final$V2[matching_row]
  }
}

BCC.ccss.org <- BCC.GRCh37
## remove missing ones
BCC.ccss.org <- BCC.ccss.org[!is.na(BCC.ccss.org$BIM_CCSS_org),]
## extract ccss_org variants
write.table(BCC.ccss.org$BIM_CCSS_org, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/extract_ccss_org_vars", col.names = F, row.names = F, quote = F)

BCC.ccss.org <- BCC.ccss.org[c("CHROM", "POS_GRCh37", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")]


write.table(BCC.ccss.org, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/all_BCC_effect_ccss_org.dat", col.names = F, row.names = F, quote = F)
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/prs_prep_data.RData")


