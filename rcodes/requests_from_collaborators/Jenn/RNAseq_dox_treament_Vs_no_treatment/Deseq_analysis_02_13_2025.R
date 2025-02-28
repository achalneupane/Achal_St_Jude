rm(list=ls())
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/eQTL_RNAseq_NW_12_6_2024/")



library(ggplot2)
library(ggfortify)
library(ggrepel)
# BiocManager::install("DESeq2")
library(DESeq2)
counts<-read.table(file='SAPKO-842478-UNSTRANDED_RSEM_gene_count.2025-01-07_19-54-28.txt', header=TRUE)
dim(counts)
# [1] 60754    28

library(dplyr)

counts$geneSymbol[duplicated(counts$geneSymbol)] <- paste0(counts$geneSymbol, "_", counts$geneID)[duplicated(counts$geneSymbol)]
# counts = dplyr::select(counts, -c(geneID, bioType, annotationLevel)) # If using gene Symbol
counts = dplyr::select(counts, -c(geneSymbol, bioType, annotationLevel)) # If using ensemble ID
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


sample_names <- colnames(counts.matrix_ROBI)
treatment <- factor(gsub(".*_(\\d+)$", "\\1", sample_names), levels = c("0", "1", "3"))  # Setting "0" as control
colData <- data.frame(row.names = sample_names, treatment = treatment)

dds <- DESeqDataSetFromMatrix(countData = counts.matrix_ROBI,
                              colData = colData,
                              design = ~ treatment)  # Model treatment effect


dds <- DESeq(dds)

res_0_vs_1 <- results(dds, contrast = c("treatment", "0", "1"))
res_0_vs_3 <- results(dds, contrast = c("treatment", "0", "3"))
res_1_vs_3 <- results(dds, contrast = c("treatment", "1", "3"))



## 1
res_0_vs_1 <- res_0_vs_1[order(res_0_vs_1$padj),]
res_0_vs_1.sub <- subset(res_0_vs_1, pvalue < 0.05)


res_0_vs_1.sub$expression = ifelse(res_0_vs_1.sub$padj < 0.05 & (res_0_vs_1.sub$log2FoldChange) >= (2), "Up",
                                                 ifelse(res_0_vs_1.sub$padj < 0.05 & (res_0_vs_1.sub$log2FoldChange) <= (-2),
                                                        "Down", "Stable"))
res_0_vs_1.sub$expression = factor(res_0_vs_1.sub$expression, levels = c("Up", "Down", "Stable"))
res_0_vs_1.sub = subset(res_0_vs_1.sub, !is.na(padj))


## 2
res_0_vs_3 <- res_0_vs_3[order(res_0_vs_3$padj),]
res_0_vs_3.sub <- as.data.frame(subset(res_0_vs_3, pvalue < 0.05))
dim(res_0_vs_3.sub)


res_0_vs_3.sub$expression = ifelse(res_0_vs_3.sub$padj < 0.05 & (res_0_vs_3.sub$log2FoldChange) >= (2), "Up",
                                   ifelse(res_0_vs_3.sub$padj < 0.05 & (res_0_vs_3.sub$log2FoldChange) <= (-2),
                                          "Down", "Stable"))
res_0_vs_3.sub$expression = factor(res_0_vs_3.sub$expression, levels = c("Up", "Down", "Stable"))
res_0_vs_3.sub = subset(res_0_vs_3.sub, !is.na(padj))

## 3
res_1_vs_3 <- res_1_vs_3[order(res_1_vs_3$padj),]
res_1_vs_3.sub <- as.data.frame(subset(res_1_vs_3, pvalue < 0.05))
dim(res_1_vs_3.sub)


res_1_vs_3.sub$expression = ifelse(res_1_vs_3.sub$padj < 0.05 & (res_1_vs_3.sub$log2FoldChange) >= (2), "Up",
                                   ifelse(res_1_vs_3.sub$padj < 0.05 & (res_1_vs_3.sub$log2FoldChange) <= (-2),
                                          "Down", "Stable"))
res_1_vs_3.sub$expression = factor(res_1_vs_3.sub$expression, levels = c("Up", "Down", "Stable"))
res_1_vs_3.sub = subset(res_1_vs_3.sub, !is.na(padj))

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/eQTL_RNAseq_NW_12_6_2024/for_Jenn_2_13_2025")

write.table(as.data.frame(res_0_vs_1.sub), "DEGs_0_vs_1.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(as.data.frame(res_0_vs_3.sub), "DEGs_0_vs_3.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(as.data.frame(res_1_vs_3.sub), "DEGs_1_vs_3.txt", sep = "\t", quote = F, row.names = T, col.names = T)



#Create Volcano Plot from DESEQ results
p = ggplot(data = res_0_vs_1.sub, aes(x = log2FoldChange,
                                  y = -log10(padj),
                                  colour = expression,
                                  label = rownames(res_0_vs_1.sub))) +
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

tiff(file="res_0_vs_1.sub.tiff",
     width = 18, height = 16, units = 'cm', res = 300, compression = 'lzw')
p

#######################################
dds_ROBIS<-DESeqDataSetFromMatrix(countData = counts.matrix_ROBI,colData = colData, design=~1)
rld_ROBIS<-vst(dds_ROBIS,blind = FALSE)
mat_ROBIS<-assay(rld_ROBIS)
pcaData<-plotPCA(rld_ROBIS, intgroup=c("treatment"), ntop=100, returnData=TRUE)
percentVar<-round(100 * attr(pcaData, "percentVar"))

library(ggplot2)
g <- ggplot(pcaData, aes(PC1, PC2))
g <- g + geom_point(size = 5, aes(col = treatment)) +
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
g

# Save the plot as a TIFF file
tiff(
  file = "Jenn_top_100_genes.tiff",
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

