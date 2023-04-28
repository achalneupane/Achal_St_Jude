########################################
# Pre-processing of raw file metadata ##
########################################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/")
md5.fastq.downloaded <- read.table("all.checked.md5", sep = "\t", header = F)
md5.fastq.downloaded <- md5.fastq.downloaded[-2]

md5.downloaded <- read.table("./md5_folder/all.md5.txt", header = F)

dim(md5.fastq.downloaded)
dim(md5.downloaded)
md5.downloaded$V2 <- gsub("./", "", md5.downloaded$V2)

# This means that all files were downloaded correctly based on MD5
sum(md5.downloaded$V2 == md5.fastq.downloaded$V3) # 120
sum(md5.downloaded$V1 == md5.fastq.downloaded$V1) # 120


## Create SRM table
path <- "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/UMD_sRNA_seq/"

R1 <- read.table("R_1_files.txt")
R1$tmpfq <- gsub("-", "_", R1$V1)
R1$ID <- sapply(strsplit(R1$tmpfq, "_"), function(x) paste(x[1:2], collapse = "_"))
R1$FASTQ1 <- paste0(path, R1$V1)
R1 <- R1[c("ID", "FASTQ1")]

R2 <- read.table("R_2_files.txt")
R2$tmpfq <- gsub("-", "_", R2$V1)
R2$ID <- sapply(strsplit(R2$tmpfq, "_"), function(x) paste(x[1:2], collapse = "_"))
R2$FASTQ2 <- paste0(path, R2$V1)
R2 <- R2[c("ID", "FASTQ2")]

sum(R1$ID != R2$ID)

df_combined <- merge(R1, R2, by = "ID")
# library(gdata)
# df_combined <- interleave(R1, R2)

df_combined$subject_ID <- sapply(strsplit(df_combined$ID,"_"), `[`, 1)

write.table(df_combined, "SRM_request.txt", col.names = T, row.names = F, sep = "\t")