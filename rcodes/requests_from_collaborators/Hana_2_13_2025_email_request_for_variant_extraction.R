# # Hana's Email on: 12/2/2024
# Hi Achal,
# Happy Monday! Hope you had a great holiday!
#   I had a question about the WGS for some of the samples we sent you before. Would it be possible to send me the list of variants for the following samples (not sure if 25-3 and 26-3 were included in the ones we sent you)?

# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_110_samples/VQSR_ApplyRecalibration_edited/HANA_Northwestern_request")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//WGS_Northwestern_Joint_call_153_samples/Hana_request_12_2_2024/")

library(dplyr)
library(data.table)
library(readxl)

# Read the first sheet
genes.hana <- read_excel("List of cardiac genes.xlsx", sheet = 1)
genes.hana$Genes <- gsub("\\(.*", "", genes.hana$Genes)
genes.hana$Genes <- gsub("\u00A0", " ", genes.hana$Genes)
genes.hana$Genes <- trimws(genes.hana$Genes)
genes.hana$Genes[grepl("TAZ", genes.hana$Genes)] <- "TAFAZZIN"

# Install if necessary
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")

# Load library
library(biomaRt)

# Connect to Ensembl (GRCh38)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 38)

# List of gene names
gene_list <- genes.hana$Genes  # Make sure this is a character vector

# Query Ensembl for start and end positions
gene_positions <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                        filters = "hgnc_symbol",
                        values = gene_list,
                        mart = ensembl)

# View results
gene_positions <- gene_positions[!grepl("PATCH|HSC", gene_positions$chromosome_name),]
gene_positions
gene_positions
gene_positions <- gene_positions[gene_positions$hgnc_symbol %in% genes.hana$Genes,]
# write.table(gene_positions, "gene_regions.txt", col.names = T, row.names = F, sep = " ", quote = F)



## Read annotated variants
annotated <- fread("combined_FIELDS2-simple.txt", header = T)
annotated <- cbind.data.frame(CHR=annotated$CHROM, POS=annotated$POS, SNP=annotated$ID, REF=annotated$REF, ALT=annotated$ALT, IMPACT=annotated$`ANN[*].IMPACT`, EFFECT=annotated$`ANN[*].EFFECT`,
                              GENE=annotated$`ANN[*].GENE`, Clinvar=annotated$CLNSIG, HGVS_C=annotated$`ANN[*].HGVS_C`, HGVS_P=annotated$`ANN[*].HGVS_P`, 
                              AF_nfe=annotated$AF_nfe, AF_afr=annotated$AF_afr, AF_eas=annotated$AF_eas, AF_sas=annotated$AF_sas, AF_fin=annotated$AF_fin, AF=annotated$AF) 

rs_numbers <- regmatches(annotated$SNP, regexpr("rs\\d+", annotated$SNP))
annotated$RsID <- sub("^(.*?;.*?);.*$", "\\1", annotated$SNP)
annotated$RsID[!grepl("rs", annotated$RsID)] <- NA
annotated$RsID <- sub("^chr[^;]+;", "", annotated$RsID)
annotated$ID <- annotated$SNP
annotated$SNP <- sub(";.*", "", annotated$SNP)

MYH7 <- annotated[grepl("MYH7", annotated$GENE),]
MYH7 <- MYH7[!duplicated(MYH7$SNP),]
# chr14:23424119:G:A

annotated <- annotated[!grepl("MYH7", annotated$GENE),]
LOF <- annotated[grepl("splice_donor_variant|splice_acceptor_variant|missense|frameshift_variant|stop_gained|inframe_deletion|inframe_insertion", annotated$EFFECT, ignore.case = T),]
clinvar <- annotated[(grepl("^Pathogenic|^Likely_pathogenic", annotated$Clinvar, ignore.case = T)),]
clinvar <- clinvar[!duplicated(clinvar$SNP),]
LOF <- LOF[!LOF$SNP %in% clinvar$SNP,]
LOF <- LOF[!duplicated(LOF$SNP),]
annotated <- rbind.data.frame(MYH7, clinvar, LOF)
annotated <- annotated[annotated$GENE %in% genes.hana$Genes,]

write.table(as.data.frame(annotated$SNP), "extract_hana_variants.txt", quote = F, row.names = F, col.names = F)

#########
## RAW ##
#########
raw <- read.delim("WGS_153_samples_hana_12_02_2024_recodeA.raw", sep = " ", header = T, check.names = F)
dim(raw)
head(raw)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

# exclude chr9.84286010.A.G
raw <- raw[!grepl("FID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]
colnames(raw) # check this with the .bim REF and NON-REFERENCE alleles

rownames(raw) <- raw$IID
raw <- raw[!grepl("IID", colnames(raw))]

df.list <- list()
for (i in 1:ncol(raw)){
  # i=1
df.list.tmp <- raw[i]

REF = unlist(strsplit(colnames(df.list.tmp), "\\:"))[3]
ALT = unlist(strsplit(colnames(df.list.tmp), "\\:"))[4]

df.list.tmp$genotype <- as.character(df.list.tmp[,1])
df.list.tmp$genotype <- gsub("0", paste(REF,REF, sep = "/"), df.list.tmp$genotype)
df.list.tmp$genotype <- gsub("1", paste(REF,ALT, sep = "/"), df.list.tmp$genotype)
df.list.tmp$genotype <- gsub("2", paste(ALT,ALT, sep = "/"), df.list.tmp$genotype)
df.list.tmp$samples <- row.names(df.list.tmp)
colnames(df.list.tmp)
df.list.tmp[1] <- df.list.tmp$genotype
df.list.tmp <- df.list.tmp[-2]
df.list[[i]] <- df.list.tmp
}

cc <- Reduce(function(...) merge(..., by= "samples", all.x=T), df.list)
dd <- t(cc)
colnames(dd) <- dd[1,]
dd <- dd[-1,]
dd <- cbind.data.frame(SNP_ID = rownames(dd), dd)

table(dd$SNP_ID %in% annotated$SNP)
dd <- cbind.data.frame(annotated, dd[match(annotated$SNP, dd$SNP_ID),])

# Flag VQSR PASS
VQSR.PASS <- fread("VQSR_PASS.txt")
dd$VQSR_PASS <- dd$SNP %in% VQSR.PASS$V2

write.table(dd, "Hana_genotype_data_12_2_2024.txt", row.names = F, col.names = T, sep = "\t", quote = F)
