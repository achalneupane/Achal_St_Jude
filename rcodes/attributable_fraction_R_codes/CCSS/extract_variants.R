
all.cancers <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", sep = "\t", header = T, stringsAsFactors = F)

# Thyroid 
df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/thyroid/thyroid.score", stringsAsFactors = F)
library(tidyr)
df <- separate(data = df, col = V1, into = c("CHROM", "POS_GRCh38", "REF", "Effect_allele"), sep = ":")
df[df$Effect_allele != df$V2, c("REF", "Effect_allele")] <- df[df$Effect_allele != df$V2, c("Effect_allele", "REF")]
df <- df[,c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "V3")]
colnames(df)[colnames(df) == "V3"] <- "Effect_size"
df$TYPE <- "PGS001354_thyroid"
df$Cancer <- "Thyroid"
df$Significant_YN <- "Y"
df$CHROM <- as.numeric(gsub("chr", "", df$CHROM))


all.cancers <- rbind.data.frame(all.cancers, df)

all.cancers.INDELS <- all.cancers[nchar(all.cancers$REF) != 1 | nchar(all.cancers$Effect_allele) != 1 ,]
nrow(all.cancers.INDELS)

library(stringr) # for str_remove function
fun <- function(a, b){
  a1 <- substr(a,1,1)
  b1 <- substr(b, 1, 1)
  d <- asplit(cbind(a, b), 1)
  ifelse(a1==b1, Recall(str_remove(a,a1), str_remove(b, b1)), d)
}

all.cancers.INDELS[c('REF', 'Effect_allele')] <-  do.call(rbind, fun(all.cancers.INDELS$REF, all.cancers.INDELS$Effect_allele))
all.cancers.INDELS


PRS.INDELS <- cbind.data.frame(CHROM = all.cancers.INDELS$CHROM, START = all.cancers.INDELS$POS_GRCh38 - (nchar(all.cancers.INDELS$REF) + nchar(all.cancers.INDELS$Effect_allele)), END = all.cancers.INDELS$POS_GRCh38)

all.cancers.SNVS <- all.cancers[!(nchar(all.cancers$REF) != 1 | nchar(all.cancers$Effect_allele) != 1),]
nrow(all.cancers.SNVS)
# 2241578

PRS.SNVS <- cbind.data.frame(CHROM = all.cancers.SNVS$CHROM, START = all.cancers.SNVS$POS_GRCh38-1, END = all.cancers.SNVS$POS_GRCh38)


PRS.INDELS$KEY <- paste(PRS.INDELS$CHROM, PRS.INDELS$START, PRS.INDELS$END, sep = ":")
PRS.INDELS <- PRS.INDELS[!duplicated(PRS.INDELS$KEY),1:3]

PRS.SNVS$KEY <- paste(PRS.SNVS$CHROM, PRS.SNVS$START, PRS.SNVS$END, sep = ":")
PRS.SNVS <- PRS.SNVS[!duplicated(PRS.SNVS$KEY),1:3]


write.table(PRS.INDELS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/Indels_PRS.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(PRS.SNVS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/SNVs_PRS.bed", row.names = F, col.names = F, quote = F, sep = "\t")
