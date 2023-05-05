##################
## CCSS_exp_WGS ##
##################

# setwd("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/Aron_Onerup_05_04_2023/")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WGS/Aron_Onerup_05_04_2023")

library(dplyr)

raw <- read.delim("aron_bfile_recodeA.raw", sep = " ", header = T, check.names = F)
dim(raw)
head(raw)
colnames(raw)[7] <- "chr5:609978:C:T"

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER


rownames(raw) <- raw$IID
raw <- raw[!grepl("IID", colnames(raw))]

df.list <- list()
for (i in 1:ncol(raw)){
  # i=1
df.list.tmp <- raw[i]

REF = unlist(strsplit(colnames(df.list.tmp), "\\:"))[3]
ALT = unlist(strsplit(colnames(df.list.tmp), "\\:"))[4]

df.list.tmp$genotype <- as.character(df.list.tmp[,1])
df.list.tmp$genotype <- gsub("0", paste0(REF,"/",REF), df.list.tmp$genotype)
df.list.tmp$genotype <- gsub("1", paste0(REF, "/", ALT), df.list.tmp$genotype)
df.list.tmp$genotype <- gsub("2", paste0(ALT, "/",ALT), df.list.tmp$genotype)
df.list.tmp$samples <- row.names(df.list.tmp)
colnames(df.list.tmp)
df.list.tmp[1] <- df.list.tmp$genotype
df.list.tmp <- df.list.tmp[-2]
df.list[[i]] <- df.list.tmp
}

cc <- Reduce(function(...) merge(..., by= "samples", all.x=T), df.list)
colnames(cc)
cc <- cc[c("samples", "chr5:609978:C:T")]



df <- read.table("QCed_samples")

sum(cc$samples %in% df$V1 )
# cc$VQSR_samples <- ifelse(cc$samples %in% df$V1, "PASS", "FAIL")

write.table(cc, "Aron_genotype_rs924607.txt", sep = "\t", col.names = T, quote = F, row.names = F)

