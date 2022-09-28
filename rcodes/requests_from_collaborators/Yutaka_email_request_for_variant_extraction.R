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

df.list <- list()
for (i in 1:ncol(df)){
  i=1
df.list <- cbind.data.frame(sapply(strsplit(df$chr12.53211761.G.A, "="), `[`, 1), sapply(strsplit(df$chr12.53211761.G.A, "="), `[`, 2))

colnames(df.list) <- c(colnames(df)[i], "genotype")
REF = unlist(strsplit(colnames(df)[i], "\\."))[3]
ALT = unlist(strsplit(colnames(df)[i], "\\."))[4]

}
