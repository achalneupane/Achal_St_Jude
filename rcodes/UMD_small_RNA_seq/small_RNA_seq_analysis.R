setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/UMD_sRNA_seq/")


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



## Clustering
data <- structure(list(K3_rep1 = c(100L, 50L, 500L, 75L, 200L), 
                       K3_rep2 = c(200L, 75L, 550L, 100L, 225L), 
                       K12_rep1 = c(150L, 100L, 600L, 50L, 250L), 
                       K12_rep2 = c(250L, 125L, 575L, 75L, 225L)), 
                  .Names = c("K3_rep1", "K3_rep2", "K12_rep1", "K12_rep2"), 
                  row.names = c("miRNA1", "miRNA2", "miRNA3", "miRNA4", "miRNA5"), 
                  class = "data.frame")

# Calculate distance matrix
dist_mat <- dist(t(data), method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(dist_mat, method = "complete")

# Plot dendrogram
plot(hc)


##################
# Example miRNA expression data
# Load the DESeq2 library
library(DESeq2)

library(DESeq2)

# Read in the count data
count_data <- structure(list(K3_rep1 = c(100L, 50L, 500L, 75L, 200L), 
                             K3_rep2 = c(200L, 75L, 550L, 100L, 225L), 
                             K12_rep1 = c(150L, 100L, 600L, 50L, 250L), 
                             K12_rep2 = c(250L, 125L, 575L, 75L, 225L)), 
                        .Names = c("K3_rep1", "K3_rep2", "K12_rep1", "K12_rep2"), 
                        row.names = c("miRNA1", "miRNA2", "miRNA3", "miRNA4", "miRNA5"), 
                        class = "data.frame")

# Create a sample information data frame
sample_info <- data.frame(sample_name = colnames(count_data),
                          condition = c(rep("K3", 2), rep("K12", 2)),
                          replicate = c(1, 2, 1, 2))

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)

# Print the top 5 differentially expressed miRNAs
head(res[order(res$padj),], 5)

#########################
## With normalization
counts <- structure(list(Sample1 = c(100, 200, 300, 400), Sample2 = c(50, 150, 250, 350), Sample3 = c(75, 175, 275, 375)), .Names = c("Sample1", "Sample2", "Sample3"), row.names = c("miR-1", "miR-2", "miR-3", "miR-4"), class = "data.frame")
coldata <- data.frame(condition = c("Control", "Control", "Treatment"))

# Calculate the library size for each sample
library_size <- colSums(counts)

# Calculate the normalization factor for each sample
median_library_size <- median(library_size)
normalization_factor <- median_library_size / library_size

# Multiply each count in the sample by its normalization factor
normalized_counts <- t(t(counts) * normalization_factor)

dds <- DESeqDataSetFromMatrix(countData = normalized_counts, colData = coldata, design = ~ condition)

# Perform differential expression analysis
dds <- DESeq(dds)
results <- results(dds)

summary(results)
