setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/Deseq_data/")

library(ggplot2)
library(ggfortify)
library(ggrepel)
library(dplyr)
library(data.table)


counts<-read.table(file='SAPKO-sjlife_128_rnaseq-STRANDED_RSEM_gene_count.2023-03-02_06-13-263_replaced_low_read_samples.txt', header=TRUE)
counts = select(counts, -c(geneSymbol, bioType, annotationLevel))
counts = counts[rowSums(counts[,-1])>40,]
counts.matrix_ROBI<-as.matrix(counts)
counts2<-counts[,-1]
rownames(counts2)<-counts[,1]
counts2<-counts2[,order(colnames(counts2))]
counts.matrix_ROBI<-as.matrix(counts2)
mode(counts.matrix_ROBI)<-"integer"

properties<-read.table('final_phenotype_file_matched_design_with_TBID_128_samples_with_RNAseq_with_anthra_jco_dose_any3_categories2.txt', sep = '\t', header=TRUE)
row.names(properties)<-properties$sjlid
properties<-properties[,-1]
properties$pair.id=as.factor(properties$pair.id)
properties$status=as.factor(properties$status)
properties<-properties[order(row.names(properties)),]
properties$anthra_jco_dose_any=as.factor(properties$anthra_jco_dose_any)
properties$agedx=as.factor(properties$agedx)
properties$age_sample=as.factor(properties$age_sample)
properties$gender=as.factor(properties$gender)
properties$racegrp=as.factor(properties$racegrp)

library(DESeq2)
dds_ROBIS<-DESeqDataSetFromMatrix(countData = counts.matrix_ROBI,colData = properties, design=~1)
rld_ROBIS<-vst(dds_ROBIS,blind = FALSE)
mat_ROBIS<-assay(rld_ROBIS)
View(properties)
pcaData<-plotPCA(rld_ROBIS, intgroup=c("status","gender"), ntop=100, returnData=TRUE)
percentVar<-round(100 * attr(pcaData, "percentVar"))
g <- ggplot(pcaData,aes(PC1,PC2))
g<-g + geom_point(size=5,aes(shape=status,col=gender,fill=gender)) + scale_shape_manual(values=21:22) + theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("Principal Component 1") + ylab("Principal Component 2") + theme(axis.text = element_text(face = "bold",color="black"),panel.border = element_blank(),axis.line = element_line(colour = "black"))+theme(axis.title.x = element_text(face = "bold"),axis.title.y = element_text(face = "bold"))+theme(axis.title.x = element_text(colour = "black"),axis.title.y = element_text(colour = "black"))+theme(text=element_text(color="black",face="bold"))
tiff(file="PCA_sjlife_128.tiff",width = 16, height = 14, units = 'cm', res = 300, compression = 'lzw')
show(g)
dev.off()


##################################################################################
# Deseq Analysis

ddset = DESeqDataSetFromMatrix(
  countData = counts.matrix_ROBI, # raw counts
  colData = properties, # phenotype data
  design = ~status+agedx+gender+racegrp+anthra_jco_dose_any+age_sample # analysis model
)
deseq_object = DESeq(ddset)
deseq_results = data.frame(results(deseq_object, contrast = c("status", "1", "0")))
deseq_results = deseq_results[order(deseq_results$padj),]
deseq_results_filtered = subset(deseq_results, padj<0.05)
dim(deseq_results_filtered)
deseq_results$expression = ifelse(deseq_results$padj < 0.05 & (deseq_results$log2FoldChange) >= (2), "Up",
                                  ifelse(deseq_results$padj < 0.05 & (deseq_results$log2FoldChange) <= (-2),
                                         "Down", "Stable"))
deseq_results$expression = factor(deseq_results$expression, levels = c("Up", "Down", "Stable"))
# deseq_results = subset(deseq_results, abs(log2FoldChange)<5 & !is.na(padj))
deseq_results = subset(deseq_results, !is.na(padj))

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
  scale_x_continuous(breaks = seq(-23, 23,3), limits = c(-23, 23))
tiff(file="Volcano_plot_sjlife_128.tiff",
     width = 18, height = 16, units = 'cm', res = 300, compression = 'lzw')
p
dev.off()

##################################################################################
# **All Black Samples**
# Joining 90 black samples to 12 new black samples and analyzing them together:
  
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(dplyr)
library(data.table)


counts<-read.table(file='SAPKO-SJLIFE-RNAseq_gene_counts_90_plus_12_new_samples.txt', header=TRUE)
counts = select(counts, -c(geneSymbol, bioType, annotationLevel))
counts = counts[rowSums(counts[,-1])>40,]
counts.matrix_ROBI<-as.matrix(counts)
counts2<-counts[,-1]
rownames(counts2)<-counts[,1]
counts2<-counts2[,order(colnames(counts2))]
counts.matrix_ROBI<-as.matrix(counts2)
mode(counts.matrix_ROBI)<-"integer"

properties<-read.table('phenotype_90_samples_plus_12_new_black_samples.txt', sep = '\t', header=TRUE)
row.names(properties)<-properties$sjlid
properties<-properties[,-1]
properties$status=as.factor(properties$status)
properties<-properties[order(row.names(properties)),]
properties$anthra_jco_dose_any=as.factor(properties$anthra_jco_dose_any)
properties$agedx=as.factor(properties$agedx)
properties$sample_age=as.factor(properties$sample_age)
properties$gender=as.factor(properties$gender)


library(DESeq2)
dds_ROBIS<-DESeqDataSetFromMatrix(countData = counts.matrix_ROBI,colData = properties, design=~1)
rld_ROBIS<-vst(dds_ROBIS,blind = FALSE)
mat_ROBIS<-assay(rld_ROBIS)
View(properties)
pcaData<-plotPCA(rld_ROBIS, intgroup=c("status","gender"), ntop=100, returnData=TRUE)
percentVar<-round(100 * attr(pcaData, "percentVar"))
g <- ggplot(pcaData,aes(PC1,PC2))
g<-g + geom_point(size=5,aes(shape=status,col=gender,fill=gender)) + scale_shape_manual(values=21:22) + theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("Principal Component 1") + ylab("Principal Component 2") + theme(axis.text = element_text(face = "bold",color="black"),panel.border = element_blank(),axis.line = element_line(colour = "black"))+theme(axis.title.x = element_text(face = "bold"),axis.title.y = element_text(face = "bold"))+theme(axis.title.x = element_text(colour = "black"),axis.title.y = element_text(colour = "black"))+theme(text=element_text(color="black",face="bold"))
tiff(file="PCA_sjlife_90_plus_12_black_samples.tiff",width = 16, height = 14, units = 'cm', res = 300, compression = 'lzw')
show(g)
dev.off()

ddset = DESeqDataSetFromMatrix(
  countData = counts.matrix_ROBI, # raw counts
  colData = properties, # phenotype data
  design = ~status+agedx+gender+anthra_jco_dose_any+sample_age # analysis model
)
'''
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
-- replacing outliers and refitting for 973 genes
-- DESeq argument 'minReplicatesForReplace' = 7 
-- original counts are preserved in counts(dds)
estimating dispersions
fitting model and testing
'''

deseq_object = DESeq(ddset)
deseq_results = data.frame(results(deseq_object, contrast = c("status", "1", "0")))
deseq_results = deseq_results[order(deseq_results$padj),]
deseq_results_filtered = subset(deseq_results, padj<0.05)
dim(deseq_results_filtered)
deseq_results$expression = ifelse(deseq_results$padj < 0.05 & (deseq_results$log2FoldChange) >= (2), "Up",
                                  ifelse(deseq_results$padj < 0.05 & (deseq_results$log2FoldChange) <= (-2),
                                         "Down", "Stable"))
deseq_results$expression = factor(deseq_results$expression, levels = c("Up", "Down", "Stable"))
# deseq_results = subset(deseq_results, abs(log2FoldChange)<5 & !is.na(padj))
deseq_results = subset(deseq_results, !is.na(padj))

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
  scale_x_continuous(breaks = seq(-31, 31,3), limits = c(-31, 31))
tiff(file="Volcano_plot_sjlife_90_plus_12_black samples.tiff",
     width = 18, height = 16, units = 'cm', res = 300, compression = 'lzw')
p
dev.off()

################################################################################
# **All White Samples**
# Joining 94 black samples to 116 new white samples and analyzing them together:

library(ggplot2)
library(ggfortify)
library(ggrepel)
library(dplyr)
library(data.table)


counts<-read.table(file='SAPKO-SJLIFE-RNAseq_gene_counts_94_plus_116_new_white_samples.txt', header=TRUE)
counts = select(counts, -c(geneSymbol, bioType, annotationLevel))
counts = counts[rowSums(counts[,-1])>40,]
counts.matrix_ROBI<-as.matrix(counts)
counts2<-counts[,-1]
rownames(counts2)<-counts[,1]
counts2<-counts2[,order(colnames(counts2))]
counts.matrix_ROBI<-as.matrix(counts2)
mode(counts.matrix_ROBI)<-"integer"

properties<-read.table('phenotype_94_plus_116_white_samples.txt', sep = '\t', header=TRUE)
row.names(properties)<-properties$sjlid
properties<-properties[,-1]
properties$status=as.factor(properties$status)
properties<-properties[order(row.names(properties)),]
properties$anthra_jco_dose_any=as.factor(properties$anthra_jco_dose_any)
properties$agedx=as.factor(properties$agedx)
properties$sample_age=as.factor(properties$sample_age)
properties$gender=as.factor(properties$gender)


library(DESeq2)
dds_ROBIS<-DESeqDataSetFromMatrix(countData = counts.matrix_ROBI,colData = properties, design=~1)
rld_ROBIS<-vst(dds_ROBIS,blind = FALSE)
mat_ROBIS<-assay(rld_ROBIS)
View(properties)
pcaData<-plotPCA(rld_ROBIS, intgroup=c("status","gender"), ntop=100, returnData=TRUE)
percentVar<-round(100 * attr(pcaData, "percentVar"))
g <- ggplot(pcaData,aes(PC1,PC2))
g<-g + geom_point(size=5,aes(shape=status,col=gender,fill=gender)) + scale_shape_manual(values=21:22) + theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("Principal Component 1") + ylab("Principal Component 2") + theme(axis.text = element_text(face = "bold",color="black"),panel.border = element_blank(),axis.line = element_line(colour = "black"))+theme(axis.title.x = element_text(face = "bold"),axis.title.y = element_text(face = "bold"))+theme(axis.title.x = element_text(colour = "black"),axis.title.y = element_text(colour = "black"))+theme(text=element_text(color="black",face="bold"))
tiff(file="PCA_sjlife_90_plus_12_black_samples.tiff",width = 16, height = 14, units = 'cm', res = 300, compression = 'lzw')
show(g)
dev.off()

ddset = DESeqDataSetFromMatrix(
  countData = counts.matrix_ROBI, # raw counts
  colData = properties, # phenotype data
  design = ~status+agedx+gender+anthra_jco_dose_any+sample_age # analysis model
)
'''
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
-- replacing outliers and refitting for 973 genes
-- DESeq argument 'minReplicatesForReplace' = 7 
-- original counts are preserved in counts(dds)
estimating dispersions
fitting model and testing
'''

deseq_object = DESeq(ddset)
deseq_results = data.frame(results(deseq_object, contrast = c("status", "1", "0")))
deseq_results = deseq_results[order(deseq_results$padj),]
deseq_results_filtered = subset(deseq_results, padj<0.05)
dim(deseq_results_filtered)
deseq_results$expression = ifelse(deseq_results$padj < 0.05 & (deseq_results$log2FoldChange) >= (2), "Up",
                                  ifelse(deseq_results$padj < 0.05 & (deseq_results$log2FoldChange) <= (-2),
                                         "Down", "Stable"))
deseq_results$expression = factor(deseq_results$expression, levels = c("Up", "Down", "Stable"))
# deseq_results = subset(deseq_results, abs(log2FoldChange)<5 & !is.na(padj))
deseq_results = subset(deseq_results, !is.na(padj))

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
  scale_x_continuous(breaks = seq(-31, 31,3), limits = c(-31, 31))
tiff(file="Volcano_plot_sjlife_90_plus_12_black samples.tiff",
     width = 18, height = 16, units = 'cm', res = 300, compression = 'lzw')
p
dev.off()