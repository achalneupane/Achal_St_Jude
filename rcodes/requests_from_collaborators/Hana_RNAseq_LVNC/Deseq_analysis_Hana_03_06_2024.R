rm(list=ls())
setwd("Z:/ResearchHome/Groups/sapkogrp/projects//RNAseq/common/LVNC_Hana_project")

#####################
## Notes from Hana ##
#####################
# I have two sets of samples, on are cardiomyocytes (CM) and the other are
# endocardial cells (EC) and would like to see the differences within each cell
# type.
#
# For each cell type I have 12 samples coming from 4 groups which have 3
# replicates  (each replicate is one hiPSC line, so they are considered
# biological replicates). I have one control group and three cardiomyopathy
# groups. These four groups are : control, LVNC (left ventricular non-compaction
# cardiomyopathy), HCM (hypertrophic cardiomyopathy) and DCM (dilated
# cardiomyopathy).
#
# Here are sample details:
#
# Control: 19-3, 25-3, 26-3
#
# Cardiomyopathy (DCM): GW10, GW53, GW168
#
# Cardiomyopathy (HCM): GW129, GW167, GW169
#
# Cardiomyopathy (LVNC): GW30, GW64, GW159
#
#
#
# The first question is to see differences between control and cardiomyopathy
# group. The second and more important question is that what is the differences
# between LVNC and other cardiomyopathy groups. Specifically, we are interested
# in the trabecular/compaction process and how the trabeculations is different
# between LVNC and two other cardiomyopathy group.
################################################################################


library(ggplot2)
library(ggfortify)
library(ggrepel)
counts<-read.table(file='SAPKO-820324-STRANDED_RSEM_gene_count.2024-02-16_18-57-58.txt', header=TRUE)
# > dim(counts)
# [1] 60754    28

library(dplyr)

counts$geneSymbol[duplicated(counts$geneSymbol)] <- paste0(counts$geneSymbol, "_", counts$geneID)[duplicated(counts$geneSymbol)]

counts = dplyr::select(counts, -c(geneID, bioType, annotationLevel))
counts = counts[rowSums(counts[,-1])>20,]

# # Remove duplictae rows
# counts <- counts[!duplicated(counts$geneSymbol),]

# Restrict yourself to canonical chromosomes and you won't run into this issue -
# as often at least. You cannot collapse counts that map to different loci to
# the same "gene" just because HGNC and ENSEMBL name things differently. ENSEMBL
# is more unique so you should ideally pick the entries you want to keep instead
# of aggregating anything. 

# counts <- counts %>% group_by(geneSymbol) %>%
# summarise(across(everything(), sum, na.rm = TRUE))

# counts.matrix_ROBI<-as.matrix(counts)
rownames(counts)<-counts[,1]
counts2<-counts[,-1]

counts.matrix_ROBI<-as.matrix(counts2)
mode(counts.matrix_ROBI)<-"integer"

properties<-read.table('phenotype_LVNC_project.txt', sep = '\t', header=TRUE)
properties$Cell_Type=as.factor(properties$Cell_Type)
properties$Group=as.factor(properties$Group)
properties$Status=as.factor(properties$Status)

library(DESeq2)
dds_ROBIS<-DESeqDataSetFromMatrix(countData = counts.matrix_ROBI,colData = properties, design=~1)
rld_ROBIS<-vst(dds_ROBIS,blind = FALSE)
mat_ROBIS<-assay(rld_ROBIS)
pcaData<-plotPCA(rld_ROBIS, intgroup=c("Group","Cell_Type"), ntop=100, returnData=TRUE)
percentVar<-round(100 * attr(pcaData, "percentVar"))

library(ggplot2)
g <- ggplot(pcaData, aes(PC1, PC2))
g <- g + geom_point(size = 5, aes(fill = Group, col = Group, shape = Cell_Type)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(face = "bold", color = "black"),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(face = "bold", color = "black"),  # Corrected line
    text = element_text(color = "black", face = "bold")
  )
# Save the plot as a TIFF file
tiff(
  file = "PCA_LVNC_Project_top_100_genes.tiff",
  width = 16,
  height = 14,
  units = 'cm',
  res = 300,
  compression = 'lzw'
)

# Print the plot
print(g)

# Close the TIFF device
dev.off()


library(pheatmap)
properties2<-properties
rownames(properties2)<-properties[,1]
properties2<-subset(properties2, select=-c(sample_ID))
myMatrix<-SummarizedExperiment::assay(rld_ROBIS)
topVarGenes<-head(order(rowVars(myMatrix),decreasing = TRUE), 500)
myMatrixforheatmap<-myMatrix[topVarGenes,]
phm<-pheatmap(myMatrixforheatmap, scale="row", cluster_rows = TRUE, cluster_cols = TRUE, annotation_col = properties2, show_colnames = TRUE, show_rownames = FALSE)
tiff(file="Heatmap_LVNC_Project_top_500_genes.tiff",width = 18, height = 32, units = 'cm', res = 300, compression = 'lzw')
phm
dev.off()



ddset = DESeqDataSetFromMatrix(
    countData = counts.matrix_ROBI, # raw counts
    colData = properties2, # phenotype data
    design = ~Status + Cell_Type # analysis model
)

deseq_object = DESeq(ddset)
deseq_results = data.frame(results(deseq_object, contrast = c("Status", "Cardiomyopathy", "Control")))
deseq_results = deseq_results[order(deseq_results$padj),]
deseq_results_filtered = subset(deseq_results, pvalue<0.05)
dim(deseq_results_filtered)
deseq_results$expression = ifelse(deseq_results$padj < 0.05 & (deseq_results$log2FoldChange) >= (2), "Up",
                                  ifelse(deseq_results$padj < 0.05 & (deseq_results$log2FoldChange) <= (-2),
                                         "Down", "Stable"))
deseq_results$expression = factor(deseq_results$expression, levels = c("Up", "Down", "Stable"))
# deseq_results = subset(deseq_results, abs(log2FoldChange)<5 & !is.na(padj))
deseq_results = subset(deseq_results, !is.na(padj))

write.table(deseq_results, file="LVNC_project_deseq_results_Cardiomyopathy_vs_Control.txt", sep="\t")

#Create Volcano Plot from DESEQ results
p = ggplot(data = deseq_results, aes(x = log2FoldChange,
                                     y = -log10(padj),
                                     colour = expression,
                                     label = rownames(deseq_results))) +
    geom_point(alpha = 0.4, size = 3.5) +
    scale_color_manual(values=c("red", "blue","grey")) +
    geom_vline(xintercept = c(-2, 2), lty=4, col="black", lwd=0.8) +
    geom_hline(yintercept = -log10(0.05), lty=4, col="black", lwd=0.8) +
    labs(x = expression(bold(log[2]~"fold-change")),
         y = expression(bold(-log[10]~"adjusted p-value"))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=16, colour = "black", face = "bold"),
          axis.text = element_text(colour = "black", face = "bold"),
          legend.title = element_blank()) +
    scale_x_continuous(breaks = seq(-10, 10,3), limits = c(-10, 10))

tiff(file="Volcano_plot_Cardiomyopathy_vs_Control.tiff",
     width = 18, height = 16, units = 'cm', res = 300, compression = 'lzw')
p
dev.off()

###############################################################################################
## 2. 
# The second and more important question is that what is the differences
# between LVNC and other cardiomyopathy groups. Specifically, we are interested
# in the trabecular/compaction process and how the trabeculations is different
# between LVNC and two other cardiomyopathy group.
###############################################################################################
## remove controls from properties 
properties3 <- properties2[properties2$Status != "Control",]
# Define LVNC group
LVNC_samples <- c("GW30_CM", "GW30_EC", "GW64_CM", "GW64_EC", "GW159_CM", "GW159_EC")

# Create a new column in colData to indicate LVNC status
properties3$LVNC_status <- ifelse(rownames(properties3) %in% LVNC_samples, "LVNC", "Other_Cardiomyopathy")

counts.matrix_ROBI.3 <- counts.matrix_ROBI[,match(rownames(properties3), colnames(counts.matrix_ROBI))]

# Run DESeq2 analysis
ddset_lvnc_vs_others <- DESeqDataSetFromMatrix(
  countData = counts.matrix_ROBI.3,
  colData = properties3,
  design = ~ LVNC_status + Cell_Type
)

deseq_object_lvnc_vs_others <- DESeq(ddset_lvnc_vs_others)
deseq_results_lvnc_vs_others <- data.frame(results(deseq_object_lvnc_vs_others, contrast = c("LVNC_status", "LVNC", "Other_Cardiomyopathy")))
deseq_results_lvnc_vs_others <- deseq_results_lvnc_vs_others[order(deseq_results_lvnc_vs_others$padj),]
deseq_results_filtered_lvnc_vs_others <- subset(deseq_results_lvnc_vs_others, pvalue < 0.05)
dim(deseq_results_filtered_lvnc_vs_others)

# Additional filtering if needed
# deseq_results_lvnc_vs_others = subset(deseq_results_lvnc_vs_others, abs(log2FoldChange) < 5 & !is.na(padj))
# deseq_results_lvnc_vs_others = subset(deseq_results_lvnc_vs_others, !is.na(padj))


deseq_results_lvnc_vs_others$expression = ifelse(deseq_results_lvnc_vs_others$padj < 0.05 & (deseq_results_lvnc_vs_others$log2FoldChange) >= (2), "Up",
                                  ifelse(deseq_results_lvnc_vs_others$padj < 0.05 & (deseq_results_lvnc_vs_others$log2FoldChange) <= (-2),
                                         "Down", "Stable"))
deseq_results_lvnc_vs_others$expression = factor(deseq_results_lvnc_vs_others$expression, levels = c("Up", "Down", "Stable"))
# deseq_results_lvnc_vs_others = subset(deseq_results_lvnc_vs_others, abs(log2FoldChange)<5 & !is.na(padj))
deseq_results_lvnc_vs_others = subset(deseq_results_lvnc_vs_others, !is.na(padj))


# Write results to a file
write.table(deseq_results_lvnc_vs_others, file = "LVNC_project_deseq_results_LVNC_vs_Other_Cardiomyopathy.txt", sep = "\t")


# baseMean log2FoldChange     lfcSE      stat       pvalue        padj
# ENSG00000189229.11  12.34759      -4.792452 0.9413042 -5.091289 3.556372e-07 0.007925732
# ENSG00000270641.1  421.42669     -12.321998 2.5507269 -4.830779 1.359999e-06 0.015154474
# ENSG00000204963.6   65.86702      -2.932972 0.6283577 -4.667678 3.046224e-06 0.019494477
# ENSG00000285756.2  288.26841      -1.521497 0.3279710 -4.639120 3.498964e-06 0.019494477
# ENSG00000110455.13 110.23448       1.051068 0.2345796  4.480644 7.441803e-06 0.031449635
# ENSG00000170689.10 128.10599      -3.216789 0.7223837 -4.453020 8.467101e-06 0.031449635

# baseMean: The average expression of the gene across all samples. It serves as a baseline for comparison.
# 
# log2FoldChange: The log2 fold change represents the change in gene expression between the two groups being compared. In this case, it's LVNC vs. Other Cardiomyopathy.
# 
# A negative log2 fold change (-4.79, for example) indicates that the gene is downregulated in LVNC compared to Other Cardiomyopathy.
# A positive log2 fold change (1.05, for example) indicates upregulation in LVNC compared to Other Cardiomyopathy.
# lfcSE: The standard error of the log2 fold change. It provides a measure of the uncertainty associated with the estimated fold change.
# 
# stat: The Wald test statistic for the significance of the log2 fold change. Larger absolute values indicate stronger evidence against the null hypothesis.
# 
# pvalue: The p-value associated with the Wald test. It indicates the probability of observing the given results if the null hypothesis (no differential expression) is true.
# 
# padj (Adjusted p-value): The p-value adjusted for multiple testing using methods like the Benjamini-Hochberg procedure. It helps control the false discovery rate.
# 
# Now, looking at the specific genes:
# 
# ENSG00000189229.11: Highly downregulated (log2FoldChange = -4.79) in LVNC compared to Other Cardiomyopathy. The result is statistically significant with a low adjusted p-value (0.0079).
# 
# ENSG00000270641.1: Highly downregulated (log2FoldChange = -12.32) in LVNC compared to Other Cardiomyopathy. Statistically significant with a low adjusted p-value (0.0152).
# 
# ENSG00000204963.6: Downregulated (log2FoldChange = -2.93) in LVNC compared to Other Cardiomyopathy. Statistically significant with a low adjusted p-value (0.0195).
# 
# ENSG00000285756.2: Downregulated (log2FoldChange = -1.52) in LVNC compared to Other Cardiomyopathy. Statistically significant with a low adjusted p-value (0.0195).
# 
# ENSG00000110455.13: Upregulated (log2FoldChange = 1.05) in LVNC compared to Other Cardiomyopathy. Statistically significant with a low adjusted p-value (0.0314).
# 
# ENSG00000170689.10: Highly downregulated (log2FoldChange = -3.22) in LVNC compared to Other Cardiomyopathy. Statistically significant with a low adjusted p-value (0.0314).

# Create Volcano Plot
p_lvnc_vs_others <- ggplot(data = deseq_results_lvnc_vs_others, aes(x = log2FoldChange,
                                                                    y = -log10(padj),
                                                                    colour = expression,
                                                                    label = rownames(deseq_results_lvnc_vs_others))) +
  geom_point(alpha = 0.4, size = 3.5) +
  scale_color_manual(values = c("red", "blue", "grey")) +
  geom_vline(xintercept = c(-2, 2), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  labs(x = expression(bold(log[2] ~ "fold-change")),
       y = expression(bold(-log[10] ~ "adjusted p-value"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 16, colour = "black", face = "bold"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = seq(-10, 10, 3), limits = c(-10, 10))

# Save the plot as a TIFF file
tiff(file = "Volcano_plot_LVNC_vs_Other_Cardiomyopathy.tiff",
     width = 18, height = 16, units = 'cm', res = 300, compression = 'lzw')
print(p_lvnc_vs_others)
dev.off()
