library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(gplots)
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


# 1. Cross-sectional Comparisons:
########################################
## a. 3 months: Poor vs Good function ##
########################################
# Define the comparison
contrast_3P <- c("three_month_status", "Poor", "Good")
contrast_12P <- c("twelve_month_status", "Poor", "Good")

# Run DESeq2 with the contrast
dds <- DESeqDataSetFromMatrix(countData = df, colData = clinical, design = ~ three_month_status + Donor_Age + Donor_Gender_M_F)
dds <- DESeq(dds)
res_3P <- results(dds, contrast = contrast_3P)


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


# unsupervised hierarchical clustering
vsd <- varianceStabilizingTransformation(dds)
# extract normalized count data
norm_counts <- assay(vsd)
# perform hierarchical clustering
dist_matrix <- dist(t(norm_counts))
hc <- hclust(dist_matrix)
plot(hc, main="Hierarchical clustering of small RNA-seq data")


# Perform unsupervised hierarchical clustering of the samples using the heatmap.2 function from the gplots package:
# scale the data
scaled_data <- t(scale(t(counts(dds)), center = TRUE, scale = TRUE))

# create a heatmap of the samples using unsupervised hierarchical clustering
heatmap.2(scaled_data,
          Colv=FALSE,
          Rowv=TRUE,
          dendrogram="row",
          trace="none",
          margins=c(10,10),
          key=TRUE,
          keysize=1.5,
          key.title="Log2 Count",
          key.xlab="",
          cexRow=0.8,
          cexCol=0.8,
          labCol=rownames(design),
          density.info="none",
          main="Unsupervised Hierarchical Clustering Heatmap")


#########################################
## b. 12 months: Poor vs Good function ##
#########################################
dds <- DESeqDataSetFromMatrix(countData = df, colData = clinical, design = ~ twelve_month_status + Donor_Age + Donor_Gender_M_F)
dds <- DESeq(dds)
res_12P <- results(dds, contrast = contrast_12P)

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


# unsupervised hierarchical clustering
vsd <- varianceStabilizingTransformation(dds)
# extract normalized count data
norm_counts <- assay(vsd)
# perform hierarchical clustering
dist_matrix <- dist(t(norm_counts))
hc <- hclust(dist_matrix)
plot(hc, main="Hierarchical clustering of small RNA-seq data")



# 2. Longitudinal comparison: 
# Longitudinal comparison: 
#####################################
## 3 months Good vs 12 months Good ##
#####################################
threeM.g.12M.g <- clinical[grepl("Good", clinical$Overall_status),]
count_matrix <- df[colnames(df) %in% (unique(threeM.g.12M.g$Subject_ID))]

# Create DESeqDataSet object using count matrix and the properties data
properties_data <- threeM.g.12M.g

# Define DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = threeM.g.12M.g, design = ~ time_point + Donor_Age + Donor_Gender_M_F)

# Run differential expression analysis
dds <- DESeq(dds)
contrast <- c("time_point", "twelve_P", "three_P")
res_good <-  results(dds, contrast = contrast)

res_good.df <- as.data.frame(res_good)
res_good.df <- res_good.df[!is.na(res_good.df$padj),]
ggplot(res_good.df, aes(x = log2FoldChange, y = -log10(padj))) +
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

#####################################
## 3 months Poor vs 12 months Poor ##
#####################################
threeM.p.12M.p <- clinical[grepl("Poor", clinical$Overall_status),]
count_matrix <- df[colnames(df) %in% (unique(threeM.p.12M.p$Subject_ID))]

# Create DESeqDataSet object using count matrix and the properties data
properties_data <- threeM.p.12M.p

# Define DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = df, colData = threeM.p.12M.p, design = ~ time_point + Donor_Age + Donor_Gender_M_F)

# Run differential expression analysis
dds <- DESeq(dds)
contrast <- c("time_point", "twelve_P", "three_P")
res_poor <-  results(dds, contrast = contrast)

res_poor.df <- as.data.frame(res_poor)
res_poor.df <- res_poor.df[!is.na(res_poor.df$padj),]
ggplot(res_poor.df, aes(x = log2FoldChange, y = -log10(padj))) +
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

################################

