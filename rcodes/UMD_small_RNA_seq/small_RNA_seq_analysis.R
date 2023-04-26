

#########################################
## Pre-processing of raw file metadata ##
#########################################
# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/")
# md5.fastq.downloaded <- read.table("all.checked.md5", sep = "\t", header = F)
# md5.fastq.downloaded <- md5.fastq.downloaded[-2]
# 
# md5.downloaded <- read.table("./md5_folder/all.md5.txt", header = F)
# 
# dim(md5.fastq.downloaded)
# dim(md5.downloaded)
# md5.downloaded$V2 <- gsub("./", "", md5.downloaded$V2)
# 
# # This means that all files were downloaded correctly based on MD5
# sum(md5.downloaded$V2 == md5.fastq.downloaded$V3) # 120
# sum(md5.downloaded$V1 == md5.fastq.downloaded$V1) # 120
# 
# 
# ## Create SRM table
# path <- "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/UMD_sRNA_seq/"
# 
# R1 <- read.table("R_1_files.txt")
# R1$tmpfq <- gsub("-", "_", R1$V1)
# R1$ID <- sapply(strsplit(R1$tmpfq, "_"), function(x) paste(x[1:2], collapse = "_"))
# R1$FASTQ1 <- paste0(path, R1$V1)
# R1 <- R1[c("ID", "FASTQ1")]
# 
# R2 <- read.table("R_2_files.txt")
# R2$tmpfq <- gsub("-", "_", R2$V1)
# R2$ID <- sapply(strsplit(R2$tmpfq, "_"), function(x) paste(x[1:2], collapse = "_"))
# R2$FASTQ2 <- paste0(path, R2$V1)
# R2 <- R2[c("ID", "FASTQ2")]
# 
# sum(R1$ID != R2$ID)
# 
# df_combined <- merge(R1, R2, by = "ID")
# # library(gdata)
# # df_combined <- interleave(R1, R2)
# 
# df_combined$subject_ID <- sapply(strsplit(df_combined$ID,"_"), `[`, 1)
# 
# write.table(df_combined, "SRM_request.txt", col.names = T, row.names = F, sep = "\t")

#############
## Analyis ##
#############
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/")
# read miRNA count files
df1 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/CAB.4062_sRNAseq/NextFlex_runs/results/sapkogrp_NextFlex_run_hsa_miRNA.count.txt", header = T)
df2 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/CAB.4062_sRNAseq/NebNext_runs/results/sapkogrp_NebNext_runs_hsa_miRNA.count.txt", header = T)

merged_df <- merge(df1, df2, by = "Tag_name", all = TRUE)
## replace NA with zero
merged_df[is.na(merged_df)] <- 0

## Cleaned most abundant tag_seq
merged_df$Most_abundant_tag_seq <- merged_df$Most_abundant_tag_seq.y
merged_df$Most_abundant_tag_seq[is.na(merged_df$Most_abundant_tag_seq)] <- merged_df$Most_abundant_tag_seq.x[is.na(merged_df$Most_abundant_tag_seq)]

# Cleaned sncRNA type
merged_df$sncRNA_type <- merged_df$sncRNA_type.y
merged_df$sncRNA_type[is.na(merged_df$sncRNA_type)] <- merged_df$sncRNA_type.x[is.na(merged_df$sncRNA_type)]

df <- merged_df[!grepl("\\.x|\\.y|sncRNA_type", colnames(merged_df))]
colnames(df) <- gsub("\\.", "_", colnames(df))

df <- df[!grepl("Most_abundant_tag_seq", colnames(df))]

## Read clinical data
clinical <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/CAD_ClinicalInformation.txt", header = T, sep = "\t", check.names = F)
colnames(clinical) <- gsub("[_]+", "_", gsub(" |\\/|\\(", "_", colnames(clinical)))
colnames(clinical) <- gsub("\\)", "", colnames(clinical))

##########################################################

##############
## Analysis ##
##############
# Cross-sectional Comparisons:
#   3 months: Poor vs Good function
#   12 months: Poor vs Good function
# Longitudinal comparison: 
#   3 months Good vs 12 months Good
#   3 months Poor vs 12 months Poor
# If possible, could you also request an unsupervised hierarchical clustering? 



df.test <- df[1:10, 1:11]
clinical.test <- clinical[grepl("R022|R044", clinical$Subject_ID),]

library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)





## 1.
# Convert count data to matrix
countData <- as.matrix(df[,grepl("_12P", colnames(df))])

# Create a column metadata data frame
colData <- clinical[,c("Subject_ID", "12_month_status")]

# Rename column names in colData to match column names in countData
colnames(colData)[1] <- "sample"

# Remove rows with missing values in colData
colData <- colData[complete.cases(colData),]



# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData, colData, ~ sample)



# Run DESeq2 analysis
dds <- DESeq(dds)

# Get differential expression results for 3 month status
res_3mo <- results(dds, contrast=c("3_month_status", "Poor", "Good"))

# Get differential expression results for 12 month status
res_12mo <- results(dds, contrast=c("12_month_status", "Poor", "Good"))

