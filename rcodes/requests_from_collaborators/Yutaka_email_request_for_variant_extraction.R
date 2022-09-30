setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data/")


# df <- read.delim("all_genotypes_from_ccss_WGS_for_Yutaka.txt", sep = "\t", header = T)

raw <- read.delim("merged.dat.yutaka_recodeA.raw", sep = " ", header = T)
dim(raw)
head(raw)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

# exclude chr9.84286010.A.G
raw <- raw[!grepl("chr9.84286010.A.G|FID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]
colnames(raw) # check this with the .bim REF and NON-REFERENCE alleles


rownames(raw) <- raw$IID
raw <- raw[!grepl("IID", colnames(raw))]

df.list <- list()
for (i in 1:ncol(raw)){
  # i=1
df.list.tmp <- raw[i]

REF = unlist(strsplit(colnames(df.list.tmp), "\\."))[3]
ALT = unlist(strsplit(colnames(df.list.tmp), "\\."))[4]

df.list.tmp$genotype <- as.character(df.list.tmp[,1])
df.list.tmp$genotype <- gsub("0", paste0(REF,REF), df.list.tmp$genotype)
df.list.tmp$genotype <- gsub("1", paste0(REF,ALT), df.list.tmp$genotype)
df.list.tmp$genotype <- gsub("2", paste0(ALT,ALT), df.list.tmp$genotype)
df.list.tmp$samples <- row.names(df.list.tmp)
colnames(df.list.tmp)
df.list.tmp <- df.list.tmp[, c("samples", "genotype")]
df.list[[i]] <- df.list.tmp
}

cc <- bind_rows(df.list, .id = "samples")
