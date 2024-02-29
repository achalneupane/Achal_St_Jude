library(ggplot2)
library(ggfortify)
library(ggrepel)
counts<-read.table(file='SAPKO-820324-STRANDED_RSEM_gene_count.2024-02-16_18-57-58.txt', header=TRUE)
library(dplyr)
counts = select(counts, -c(geneSymbol, bioType, annotationLevel))
counts = counts[rowSums(counts[,-1])>20,]
counts.matrix_ROBI<-as.matrix(counts)
counts2<-counts[,-1]
rownames(counts2)<-counts[,1]
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
g <- ggplot(pcaData,aes(PC1,PC2))
g<-g + geom_point(size=5,aes(fill=Group, col=Group, shape=Cell_Type)) + theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("Principal Component 1") + ylab("Principal Component 2") + theme(axis.text = element_text(face = "bold",color="black"),panel.border = element_blank(),axis.line = element_line(colour = "black"))+theme(axis.title.x = element_text(face = "bold"),axis.title.y = element_text(face = "bold"))+theme(axis.title.x = element_text(colour = "black"),axis.title.y = element_text(colour = "black"))+theme(text=element_text(color="black",face="bold"))
tiff(file="PCA_LVNC_Project_top_100_genes.tiff",width = 16, height = 14, units = 'cm', res = 300, compression = 'lzw')
show(g)
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
