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

## Wanted columns
ST6 <- ST6[,c("CHR", "hg38BP", "OA", "EA", "beta")]
ST6$TYPE <- "ST6"
ST6$TYPE2 <- "ST6"
colnames(ST6) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")

###############
## PGS003416 ##
###############

PGS003416 <- read.table("PGS003416_hmPOS_GRCh38.txt", sep = "\t", header = T)
PGS003416 <- PGS003416[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight")]
PGS003416$TYPE <- "PGS003416"
PGS003416$TYPE2 <- "PGS003416"
colnames(PGS003416) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")

###############
## PGS000454 ##
###############

PGS000454 <- read.table("PGS000454_hmPOS_GRCh38.txt", sep = "\t", header = T)
PGS000454 <- PGS000454[,c("chr_name", "hm_pos", "other_allele", "effect_allele", "effect_weight")]
PGS000454$TYPE <- "PGS000454"
PGS000454$TYPE2 <- "PGS000454"
colnames(PGS000454) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size", "TYPE", "TYPE2")

###############
## PGS000356 ##
###############

PGS000356 <- read.table("PGS000356_hmPOS_GRCh38.txt", sep = "\t", header = T)
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
################################
## extract variants; create bed file
BED <- cbind.data.frame(BCC$CHROM, BCC$POS_GRCh38-2, BCC$POS_GRCh38+2)
write.table(BED, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/all_bed_BCC.bed", col.names = F, row.names = F, quote = F)

BCC$KEY <- paste0(BCC$CHROM, ":", BCC$POS_GRCh38, ":")

BIM.ALL <- {}
for (CHR in 1:22){
print(paste0("Doing CHR: ", CHR))
BIM <- read.table(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr",CHR, ".PASS.decomposed.bim"), header = F)
BIM$KEY <- paste0("chr", BIM$V1, ":", BIM$V4, ":")
BIM <- BIM[BIM$KEY %in% BCC$KEY,]
BIM.ALL <- rbind.data.frame(BIM.ALL,BIM)
}

# cc <- BCC[BCC$CHROM == "chr1",]
# cc <- cc [!duplicated(cc$KEY),]
