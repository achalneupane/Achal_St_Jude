library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)

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

colnames(clinical) <- gsub('^3', 'three', colnames(clinical))
colnames(clinical) <- gsub('^12', 'twelve', colnames(clinical))

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




clinical_3P <- clinical %>% mutate(Subject_ID = paste0(Subject_ID, "_3P"))
clinical_12P <- clinical %>% mutate(Subject_ID = paste0(Subject_ID, "_12P"))

clinical_3P$time_point <- "three_P"
clinical_12P$time_point <- "twelve_P"

clinical_alt <- rbind(clinical_3P[order(clinical_3P$Subject_ID), ],
                      clinical_12P[order(clinical_12P$Subject_ID), ])
clinical <- rbind(clinical_3P, clinical_12P)[order(c(seq(1,nrow(clinical_3P)*2,2),seq(2,nrow(clinical_12P)*2,2))),]

rownames(df) <- df$Tag_name
df <- select(df, -c(Tag_name))

df.test <- df[1:10, c(grep("R022|R044|R046|R058|R085", colnames(df)))]
clinical.test <- clinical[grepl("R022|R044|R046|R058|R085", clinical$Subject_ID), c("Subject_ID", "twelve_month_status", "Overall_status", "time_point", "Donor_Age", "Donor_Gender_M_F")]




# 1. Cross-sectional Comparisons:
# a. 3 months: Poor vs Good function
# Define the comparison
contrast_3P <- c("three_month_status", "Poor", "Good")
contrast_12P <- c("twelve_month_status", "Poor", "Good")

# Run DESeq2 with the contrast
dds <- DESeqDataSetFromMatrix(countData = df, colData = clinical, design = ~ three_month_status + Donor_Age + Donor_Gender_M_F)
dds <- DESeq(dds)
res_3P <- results(dds, contrast = contrast_3P)
# b. 12 months: Poor vs Good function
dds <- DESeqDataSetFromMatrix(countData = df, colData = clinical, design = ~ twelve_month_status + Donor_Age + Donor_Gender_M_F)
dds <- DESeq(dds)
res_12P <- results(dds, contrast = contrast_12P)


# Create the volcano plots
# Volcano plot for 3 months
res_3P.df <- as.data.frame(res_3P)
res_3P.df <- res_3P.df[!is.na(res_3P.df$padj),]
ggplot(res_3P.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "significant", "not_significant")), alpha = 5, size = 5) +
  scale_color_manual(values = c("red", "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  # geom_text(aes(label = ifelse(padj < 0.05, rownames(res_3P.df), "")), vjust = 2.5, size = 3) +
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +
  ggtitle("Volcano Plot for Small RNA seq Data (3 months)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()
  )



# Volcano plot for 12 months
res_12P.df <- as.data.frame(res_12P)
res_12P.df <- res_12P.df[!is.na(res_12P.df$padj),]
ggplot(res_12P.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "significant", "not_significant")), alpha = 5, size = 5) +
  scale_color_manual(values = c("red", "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  # geom_text(aes(label = ifelse(padj < 0.05, rownames(res_3P.df), "")), vjust = 2.5, size = 3) +
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +
  ggtitle("Volcano Plot for Small RNA seq Data (12 months)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()
  )


# 2. Longitudinal comparison: 
# Longitudinal comparison: 
#   3 months Good vs 12 months Good
threeM.g.12M.g <- clinical.test[grepl("Good", clinical.test$Overall_status),]

#   3 months Poor vs 12 months Poor




