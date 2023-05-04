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


clinical_3P <- clinical %>% mutate(Subject_ID = paste0(Subject_ID, "_3P"))
clinical_12P <- clinical %>% mutate(Subject_ID = paste0(Subject_ID, "_12P"))

clinical_3P$time_point <- "three_P"
clinical_12P$time_point <- "twelve_P"

clinical_alt <- rbind(clinical_3P[order(clinical_3P$Subject_ID), ],
                      clinical_12P[order(clinical_12P$Subject_ID), ])
clinical <- rbind(clinical_3P, clinical_12P)[order(c(seq(1,nrow(clinical_3P)*2,2),seq(2,nrow(clinical_12P)*2,2))),]

rownames(df) <- df$Tag_name
df <- select(df, -c(Tag_name))
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

# # clinical.test <- clinical[grepl("R022_|R044_|R046_|R058_|R085_", clinical$Subject_ID), c("Subject_ID", "three_month_status", "twelve_month_status", "time_point", "Donor_Age", "Donor_Gender_M_F")]
# clinical.test <- clinical[c(22,3,23,10,29,28,7,26,9,21,25,14,13,12,16,30,27,17,4,6), c("Subject_ID", "three_month_status", "twelve_month_status", "time_point", "Donor_Age", "Donor_Gender_M_F")]
# # df.test <- df[1:10, colnames(df) %in% clinical.test$Subject_ID]
# df.test <- df[1:30, colnames(df) %in% clinical.test$Subject_ID]

# 1. Cross-sectional Comparisons:
########################################
## a. 3 months: Poor vs Good function ##
########################################
# Define the comparison
contrast_3P <- c("three_month_status", "Poor", "Good")
df_filtered = df[rowSums(df)>40,]

# Run DESeq2 with the contrast
dds <- DESeqDataSetFromMatrix(countData = df_filtered, colData = clinical, design = ~ three_month_status + Donor_Age + Donor_Gender_M_F)
dds <- DESeq(dds)
res_3P <- results(dds, contrast = contrast_3P)


# Create the volcano plots
# Volcano plot for 3 months
res_3P.df <- as.data.frame(res_3P)
res_3P.df <- res_3P.df[!is.na(res_3P.df$padj),]
top20 <- res_3P.df[head(order(res_3P.df$padj, decreasing = FALSE ), 20),]
ggplot(res_3P.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "significant", "not_significant")), alpha = 5, size = 5) +
  scale_color_manual(values = c("red", "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  # geom_text(aes(label = ifelse(padj < 0.05, rownames(res_3P.df), "")), vjust = 2.5, size = 3) +
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +
  ggtitle("3 months: poor vs. good") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()
  )




#########################################
## b. 12 months: Poor vs Good function ##
#########################################
df_filtered = df[rowSums(df)>40,]
dds <- DESeqDataSetFromMatrix(countData = df_filtered, colData = clinical, design = ~ twelve_month_status + Donor_Age + Donor_Gender_M_F)
dds <- DESeq(dds)
contrast_12P <- c("twelve_month_status", "Poor", "Good")
res_12P <- results(dds, contrast = contrast_12P)

# Volcano plot for 12 months
res_12P.df <- as.data.frame(res_12P)
res_12P.df <- res_12P.df[!is.na(res_12P.df$padj),]
top20 <- res_12P.df[head(order(res_12P.df$padj, decreasing = FALSE ), 20),]
ggplot(res_12P.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "significant", "not_significant")), alpha = 5, size = 5) +
  scale_color_manual(values = c("red", "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  # geom_text(aes(label = ifelse(padj < 0.05, rownames(res_3P.df), "")), vjust = 2.5, size = 3) +
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +
  ggtitle("12 months: poor vs. good") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()
  )



# 2. Longitudinal comparison: 
# Longitudinal comparison: 
#####################################
## 3 months Good vs 12 months Good ##
#####################################
three_months_good_samples <- clinical$Subject_ID[clinical$time_point == "three_P" & clinical$three_month_status == "Good"]
twelve_months_good_samples <- clinical$Subject_ID[clinical$time_point == "twelve_P" & clinical$twelve_month_status == "Good"]
good_samples <- c(three_months_good_samples, twelve_months_good_samples)
df_filtered <- df[, good_samples]
df_filtered = df_filtered[rowSums(df_filtered)>40,]
clinical_filtered <- clinical[clinical$Subject_ID %in% good_samples,]


# Create DESeqDataSet object using count matrix and the properties data
dds <- DESeqDataSetFromMatrix(countData = df_filtered, colData = clinical_filtered, design = ~ time_point + Donor_Age + Donor_Gender_M_F)

# Run differential expression analysis
dds <- DESeq(dds)
contrast <- c("time_point", "three_P", "twelve_P")
res_good <-  results(dds, contrast = contrast)

res_good.df <- as.data.frame(res_good)
res_good.df <- res_good.df[!is.na(res_good.df$padj),]
top20 <- res_good.df[head(order(res_good.df$padj, decreasing = FALSE ), 20),]
ggplot(res_good.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "significant", "not_significant")), alpha = 5, size = 5) +
  scale_color_manual(values = c("red", "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  # geom_text(aes(label = ifelse(padj < 0.05, rownames(res_3P.df), "")), vjust = 2.5, size = 3) +
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +
  ggtitle("Good function: 3 months vs. 12 months") +
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
three_months_poor_samples <- clinical$Subject_ID[clinical$time_point == "three_P" & clinical$three_month_status == "Poor"]
twelve_months_poor_samples <- clinical$Subject_ID[clinical$time_point == "twelve_P" & clinical$twelve_month_status == "Poor"]
poor_samples <- c(three_months_poor_samples, twelve_months_poor_samples)
df_filtered <- df[, poor_samples]
df_filtered = df_filtered[rowSums(df_filtered)>40,]
clinical_filtered <- clinical[clinical$Subject_ID %in% poor_samples,]

# Create DESeqDataSet object using count matrix and the properties data
dds <- DESeqDataSetFromMatrix(countData = df_filtered, colData = clinical_filtered, design = ~ time_point + Donor_Age + Donor_Gender_M_F)

# Run differential expression analysis
dds <- DESeq(dds)
contrast <- c("time_point", "three_P", "twelve_P")
res_poor <-  results(dds, contrast = contrast)

res_poor.df <- as.data.frame(res_poor)
res_poor.df <- res_poor.df[!is.na(res_poor.df$padj),]
top20 <- res_poor.df[head(order(res_poor.df$padj, decreasing = FALSE ), 20),]
ggplot(res_poor.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "significant", "not_significant")), alpha = 5, size = 5) +
  scale_color_manual(values = c("red", "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  # geom_text(aes(label = ifelse(padj < 0.05, rownames(res_3P.df), "")), vjust = 2.5, size = 3) +
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +
  ggtitle("Poor function: 3 months vs. 12 months") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()
  )

############################################# # 4. Unsupervised hierarchical
#clustering ## ############################################ we will use no
#design parameter as specified "~1". It is a formula notation in R to indicate
#the intercept-only model, which means that there is no specific variable used
#for normalization or comparison between groups.
dds <- DESeqDataSetFromMatrix(countData = df, colData = clinical, design=~1)
dds <- DESeq(dds)

# unsupervised hierarchical clustering
vsd <- varianceStabilizingTransformation(dds)
# extract normalized count data
norm_counts <- assay(vsd)
# Perform unsupervised hierarchical clustering of the samples using the heatmap.2 function from the gplots package:
# scale the data
# Filter dds to include only the top 20 genes by log2 fold change
# top_20_genes <- head(order(results(dds)$log2FoldChange, decreasing = TRUE), 20)
# dds_top_20 <- dds[top_20_genes,]
top_genes <- rownames(norm_counts)[order(rowVars(norm_counts), decreasing = TRUE)[1:100]]
top_norm_counts <- norm_counts[top_genes, ]
# Scale the gene expression data
# scaled_data <- t(scale(t(counts(dds_top_20)), center = TRUE, scale = TRUE))

## create a heatmap of the samples using unsupervised hierarchical clustering

# Create a vector of colors based on twelve_month_status
# color_vector <- ifelse(clinical$three_month_status == "Good", "green", "red")
color_vector <- ifelse(clinical$twelve_month_status == "Good", "green", "red")



# Plot the heatmap with colored column names
heatmap.2(top_norm_counts,
          Colv=TRUE,
          Rowv=TRUE,
          distfun = dist,
          hclustfun = hclust,
          trace="none",
          margins=c(10,10),
          key=TRUE,
          keysize=1.5,
          key.title="Log2 Count",
          key.xlab="",
          cexRow=0.8,
          cexCol=0.8,
          reorderfun = function(d, w) reorder(d, w),
          density.info="none",
          main="",
          ColSideColors=color_vector)

# Define the colors for the legend
colors <- c("Good" = "green", "Poor" = "red")

# Add the legend
legend("top", legend = levels(factor(clinical$twelve_month_status)), col = colors, pch = 15, horiz = T)



