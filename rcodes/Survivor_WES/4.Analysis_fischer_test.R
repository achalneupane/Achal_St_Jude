library(data.table)
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife/phenotype.RData")
## 1. clinvar
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar <- clinvar[!duplicated(clinvar$SNP),]
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)


# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar.all <- clinvar[which(clinvar$AF <0.01),]
clinvar.eur <- clinvar[which(clinvar$AF <0.01 & clinvar$AF_nfe < 0.01),]
clinvar.afr <- clinvar[which(clinvar$AF <0.01 & clinvar$AF_afr < 0.01),]


loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)
# - 
# 70494 

loftee <- loftee[!duplicated(loftee$SNP),]


loftee$AF <- as.numeric(loftee$AF.1)
loftee$AF_nfe <- as.numeric(loftee$AF_nfe)
loftee$AF_afr <- as.numeric(loftee$AF_afr)
loftee$CHROM  <- sub("([0-9XY]+):.+", "\\1", loftee$SNP)
loftee$POS <- sub("chr[0-9XY]+:(\\d+):.+", "\\1", loftee$SNP)

# make rare
loftee.all <- loftee[which(loftee$AF <0.01),]
loftee.eur <- loftee[which(loftee$AF <0.01 & loftee$AF_nfe < 0.01),]
loftee.afr <- loftee[which(loftee$AF <0.01 & loftee$AF_afr < 0.01),]


snpeff <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/missense_variants_with_overlap_in_more_than_90_percent_of_prediction_tools_all_cols.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(snpeff)
snpeff$SNP <- sub(";.*", "", snpeff$ID)
snpeff <- snpeff[!duplicated(snpeff$SNP),]

snpeff$AF <- as.numeric(snpeff$AF)
snpeff$AF_nfe <- as.numeric(snpeff$AF_nfe)
snpeff$AF_afr <- as.numeric(snpeff$AF_afr)


# make rare
snpeff.all <- snpeff[which(snpeff$AF <0.01),]
snpeff.eur <- snpeff[which(snpeff$AF <0.01 & snpeff$AF_nfe < 0.01),]
snpeff.afr <- snpeff[which(snpeff$AF <0.01 & snpeff$AF_afr < 0.01),]


cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
# 46076

cc <- c(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))

bim.QC.sjlife.PLP <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant//sjlife/all_rare_variants_maf0.01_all_sjlife.bim")
table(cc %in% bim.QC.sjlife.PLP$V2)
# FALSE  TRUE 
# 16818 29254

raw <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant//sjlife/all_rare_variants_maf0.01_all_sjlife_recodeA.raw")
sum(duplicated(raw$IID))
dim(raw)
# 0
# raw <- raw[!(duplicated(raw$IID)),]
raw <- as.data.frame(raw)
rownames(raw) <- raw$IID
raw <- raw[,-c(1:6)]
# HEADER <- colnames(raw)[-c(1:6)]
HEADER=colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(bim.QC.sjlife.PLP$V2 %in% colnames(raw))
# TRUE 
# 29254 
# colnames(raw)[!colnames(raw) %in% bim.QC.sjlife.PLP$V2]


## run fisher exact test with grade 2 or higher ones with clinvar.all
raw.clinvar.all <- raw[colnames(raw) %in% clinvar.all$SNP ]
# dim(raw.clinvar.all)
# [1] 4516 6451

raw.clinvar.all$carrier <- ifelse(rowSums(raw.clinvar.all[grepl("chr", colnames(raw.clinvar.all))], na.rm = T) > 0, 1, 0)
table(raw.clinvar.all$carrier)
# 0    1 
# 701 3815 

CTCAE.data.4$sjlid %in% raw.clinvar.all$