rm(list=ls())
## read WES annotation files for POI option 1
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/")

## read Ke et al genes
Ke_genes <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al/Ke_et_al_genes.txt", header= T)
# Shekari_genes <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Shekari_et_al/Shekari_et_al_genes.txt", header= T, sep = "\t")
# Shekari_genes[!Shekari_genes$Gene %in% Ke_genes$Gene,]


## 1. clinvar
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)

tt <- cbind.data.frame(clinvar$SNP, clinvar$ID, clinvar$AF_nfe, clinvar$AF_afr, clinvar$AF)

# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar$new_GENE.clinvar <- gsub("-.*", "", clinvar$ANN....GENE)
clinvar <- cbind.data.frame(CHROM=clinvar$CHROM, POS=clinvar$POS, REF=clinvar$REF, ALT=clinvar$ALT, SNP=clinvar$SNP, ID=clinvar$ID, GENE=clinvar$ANN....GENE, new_GENE.clinvar=clinvar$new_GENE.clinvar, AF_nfe=clinvar$AF_nfe, AF_afr=clinvar$AF_afr, AF=clinvar$AF, Effect=clinvar$ANN....EFFECT, CLNSIG=clinvar$CLNSIG)
clinvar.save <- clinvar

Ke_genes$Gene[!Ke_genes$Gene %in% clinvar$new_GENE.clinvar]
# [1] "ANKRD31"   "ATG9A"     "BMP15"     "BMPR1A"    "BMPR1B"    "BNC1"      "C14orf39"  "DIAPH2"    "DMC1"      "EIF4ENIF1" "ESR2"      "EXO1"      "FANCL"    
# [14] "FIGLA"     "FOXE1"     "KHDRBS1"   "MCM9"      "NANOS3"    "NOG"       "NOTCH2"    "PGRMC1"    "POF1B"     "PRDM9"     "PSMC3IP"   "RAD51"     "SALL4"    
# [27] "SGO2"      "SPIDR"     "SYCE1"     "SYCP2L"    "TP63"      "TWNK"      "XRCC2"


## Keep only those in Ke et al
clinvar <- clinvar[clinvar$new_GENE.clinvar %in% Ke_genes$Gene,]

clinvar.all <- clinvar[which(clinvar$AF <0.01),]
clinvar.eur <- clinvar[which(clinvar$AF_nfe < 0.01),]
clinvar.afr <- clinvar[which(clinvar$AF_afr < 0.01),]




loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)
# - 
# 72368  


# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
loftee$new_GENE.loftee <- gsub("-.*", "", loftee$SYMBOL)

Ke_genes$Gene[!Ke_genes$Gene %in% loftee$new_GENE.loftee]
## Keep only those in ACMG
loftee <- loftee[loftee$new_GENE.loftee %in% Ke_genes$Gene,]


loftee$AF <- as.numeric(loftee$AF.1)
loftee$AF_nfe <- as.numeric(loftee$AF_nfe)
loftee$AF_afr <- as.numeric(loftee$AF_afr)
loftee$CHROM  <- sub("([0-9XY]+):.+", "\\1", loftee$SNP)
loftee$POS <- sub("chr[0-9XY]+:(\\d+):.+", "\\1", loftee$SNP)

# make rare
loftee.all <- loftee[which(loftee$AF <0.01),]
loftee.eur <- loftee[which(loftee$AF_nfe < 0.01),]
loftee.afr <- loftee[which(loftee$AF_afr < 0.01),]


snpeff <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/missense_variants_with_overlap_in_more_than_90_percent_of_prediction_tools_all_cols.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(snpeff)
snpeff$SNP <- sub(";.*", "", snpeff$ID)

snpeff$AF <- as.numeric(snpeff$AF)
snpeff$AF_nfe <- as.numeric(snpeff$AF_nfe)
snpeff$AF_afr <- as.numeric(snpeff$AF_afr)

snpeff$new_GENE.snpeff <- gsub("-.*", "", snpeff$ANN....GENE)
snpeff <- cbind.data.frame(CHROM=snpeff$CHROM, POS=snpeff$POS, REF=snpeff$REF, ALT=snpeff$ALT, SNP=snpeff$SNP, ID=snpeff$ID, GENE=snpeff$ANN....GENE, new_GENE.snpeff=snpeff$new_GENE.snpeff, AF_nfe=snpeff$AF_nfe, AF_afr=snpeff$AF_afr, AF=snpeff$AF, Effect=snpeff$ANN....EFFECT, CLNSIG=snpeff$CLNSIG)
snpeff.save <- snpeff

Ke_genes$Gene[!Ke_genes$Gene %in% snpeff$new_GENE.snpeff]
## Keep only those in Ke et al
snpeff <- snpeff[snpeff$new_GENE.snpeff %in% Ke_genes$Gene,]

# make rare
snpeff.all <- snpeff[which(snpeff$AF <0.01),]
snpeff.eur <- snpeff[which(snpeff$AF_nfe < 0.01),]
snpeff.afr <- snpeff[which(snpeff$AF_afr < 0.01),]


length(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
# 706
cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
# 706
colnames(cc) <- "SNP"

## QCed Bim
library(data.table)
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated.bim")
cc$SNP[!cc$SNP %in% bim.QC$V2]
cc.final <- cc[cc$SNP %in% bim.QC$V2,]

write.table(cc.final, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al_POI_rare_variants_to_extract.txt", row.names = F, col.names = F, quote = F)
length(cc.final)
# 580

# ## Autosomal recessive gene variants
# AR.genes <- ACMG$GENE[ACMG$inheritence == "AR"]
# AR.genes.clinvar <- clinvar$SNP[clinvar$new_GENE.clinvar %in% AR.genes]
# AR.genes.loftee <- loftee$SNP[loftee$new_GENE.loftee %in% AR.genes]
# AR.genes.snpeff <- snpeff$SNP[snpeff$new_GENE.snpeff %in% AR.genes]
# AR.variants <- (unique(c(AR.genes.clinvar, AR.genes.loftee, AR.genes.snpeff)))

## preQC
bim.preQC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim")
cc <- as.character(cc$SNP)
cc[!cc %in% bim.preQC$V2]
table(cc %in% bim.preQC$V2)
# TRUE 
# 706

## QCed for GQ, DP and VQSR only
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim")
table(cc %in% bim.QC$V2)
# FALSE  TRUE 
# 117   589 

## Final SJLIFE QCed data
table(cc %in% bim.QC$V2)
# FALSE  TRUE 
# 126   580


## Create carrier status
raw <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/POI_genetics//Ke_et_al_POI_rare_variants_ALL_recodeA.raw")
raw <- as.data.frame(raw)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
# HEADER <- colnames(raw)[-c(1:6)]
HEADER=colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(colnames(raw) %in% cc)
# TRUE 
# 580
## Looks good!

## Make AR gene variants either 0 or 2
# AR.variants.raw <- colnames(raw)[colnames(raw) %in% AR.variants]
# raw[, AR.variants.raw] <- ifelse(raw[, AR.variants.raw] == 2, 2, 0)


# save.image("Ke_et_al_POI_rare_variants_data.RData")
# load("Ke_et_al_POI_rare_variants_data.RData")


# Some ccss/sjlife overlapping samples in WGS data were renames as SJLIFE, but not in WES, so naming them back to CCSS.
sjlife_ccss_exp_overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sjlife_ccss_exp_overlaps.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/pca/Survivor_WGS_EUR_based_on_1kGP_Phase_3_data.txt", header = T)
# Find matching indices for SJLID in EUR$IID and replace with CCSSID
matches <- match(EUR$IID, sjlife_ccss_exp_overlaps$SJLID)
EUR$IID[!is.na(matches)] <- sjlife_ccss_exp_overlaps$CCSSID[matches[!is.na(matches)]]

AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/pca/Survivor_WGS_AFR_based_on_1kGP_Phase_3_data.txt", header = T)
# Find matching indices for SJLID in EUR$IID and replace with CCSSID
matches <- match(AFR$IID, sjlife_ccss_exp_overlaps$SJLID)
AFR$IID[!is.na(matches)] <- sjlife_ccss_exp_overlaps$CCSSID[matches[!is.na(matches)]]


# Read the tab-separated file
filtered_wes_samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/filtered_wes_samples_by_cohort.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
WES.fam <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated_unique.fam", header = F)
table(WES.fam$V2 %in% filtered_wes_samples$SampleID)
# FALSE  TRUE
# 272  7758

EUR$cohort <- filtered_wes_samples$Source[match(EUR$IID, filtered_wes_samples$SampleID)]
EUR <- EUR[!grepl("African", EUR$cohort),]
EUR <- EUR[!is.na(EUR$cohort),]
table(EUR$cohort)
# ccss_exp community.controls             sjlife
# 2299                375               3345

AFR$cohort <- filtered_wes_samples$Source[match(AFR$IID, filtered_wes_samples$SampleID)]
AFR <- AFR[!is.na(AFR$cohort),]
table(AFR$cohort)
# ccss_exp community.controls             sjlife   Survivor_African
# 100                 31                642                 76

## EUR.ccss <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSSEXP_EUR_top_20_PCs.eigenvec.ccssid", header = F)
## AFR.ccss <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSSEXP_AFR_top_20_PCs.eigenvec.ccssid", header = F)
## EUR.sjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/final_EUR-PCAS_EUR.eigenvec", header = F)
## AFR.sjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/final_AFR-PCAS_AFR.eigenvec", header = F)

# table(EUR.sjlife$V2 %in% EUR$IID[EUR$cohort=="sjlife"])
# table(AFR.sjlife$V2 %in% AFR$IID[AFR$cohort=="sjlife"])

table(rownames(raw) %in% EUR$IID)
# FALSE  TRUE
# 2019  6001

## Extract European and African genotype data
raw.eur <- raw[rownames(raw) %in% EUR$IID,]
raw.afr <- raw[rownames(raw) %in% AFR$IID,]
dim(raw.eur)
# [1] 6001  728
dim(raw.afr)
# [1] 846 728
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/ACMG/utils.R")



##################################
## 1.a Extract clinvar European ##
##################################
raw.clinvar.eur <- raw.eur[which(colnames(raw.eur) %in% clinvar.eur$SNP)]
clinvar.eur <- clinvar.eur[clinvar.eur$SNP %in% colnames(raw.eur),]

genes <- unique(clinvar.eur$new_GENE.clinvar)
carriers.clinvar.eur <- setNames(as.data.frame(rownames(raw.eur)), "IID")
snps.with.carriers.clinvar.eur <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(clinvar.eur$SNP[clinvar.eur$new_GENE.clinvar == wanted.gene])
  raw.wanted <- raw.clinvar.eur[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.clinvar.eur <- c(snps.with.carriers.clinvar.eur, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.clinvar.eur[wanted.gene] <- carriers.tmp$carrier[match(carriers.clinvar.eur$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.clinvar.eur)
# 256
dim(carriers.clinvar.eur)
# [1] 6001   56
## Adding carrier status for "All_Genes"
## Check if any gene column in carriers.clinvar.eur has carrier status (1) across all genes
# carriers.clinvar.eur$All_Genes <- ifelse(rowSums(carriers.clinvar.eur[ , -1], na.rm = TRUE) > 0, 1, 0)
dim(carriers.clinvar.eur)

#################################
## 1.b Extract clinvar African ##
#################################
raw.clinvar.afr <- raw.afr[which(colnames(raw.afr) %in% clinvar.afr$SNP)]
clinvar.afr <- clinvar.afr[clinvar.afr$SNP %in% colnames(raw.afr),]

genes <- unique(clinvar.afr$new_GENE.clinvar)
carriers.clinvar.afr <- setNames(as.data.frame(rownames(raw.afr)), "IID")
snps.with.carriers.clinvar.afr <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(clinvar.afr$SNP[clinvar.afr$new_GENE.clinvar == wanted.gene])
  raw.wanted <- raw.clinvar.afr[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.clinvar.afr <- c(snps.with.carriers.clinvar.afr, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.clinvar.afr[wanted.gene] <- carriers.tmp$carrier[match(carriers.clinvar.afr$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.clinvar.afr)
## 35
# carriers.clinvar.afr$All_Genes <- ifelse(rowSums(carriers.clinvar.afr[ , -1], na.rm = TRUE) > 0, 1, 0)

#################################
## 2.a Extract loftee European ##
#################################
raw.loftee.eur <- raw.eur[which(colnames(raw.eur) %in% loftee.eur$SNP)]
loftee.eur <- loftee.eur[loftee.eur$SNP %in% colnames(raw.eur),]

genes <- unique(loftee.eur$new_GENE.loftee)
carriers.loftee.eur <- setNames(as.data.frame(rownames(raw.eur)), "IID")
snps.with.carriers.loftee.eur <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(loftee.eur$SNP[loftee.eur$new_GENE.loftee == wanted.gene])
  raw.wanted <- raw.loftee.eur[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.loftee.eur <- c(snps.with.carriers.loftee.eur, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.loftee.eur[wanted.gene] <- carriers.tmp$carrier[match(carriers.loftee.eur$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.loftee.eur)
## 66
# carriers.loftee.eur$All_Genes <- ifelse(rowSums(carriers.loftee.eur[ , -1], na.rm = TRUE) > 0, 1, 0)

################################
## 2.b Extract loftee African ##
################################
raw.loftee.afr <- raw.afr[which(colnames(raw.afr) %in% loftee.afr$SNP)]
loftee.afr <- loftee.afr[loftee.afr$SNP %in% colnames(raw.afr),]

genes <- unique(loftee.afr$new_GENE.loftee)
carriers.loftee.afr <- setNames(as.data.frame(rownames(raw.afr)), "IID")
snps.with.carriers.loftee.afr <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(loftee.afr$SNP[loftee.afr$new_GENE.loftee == wanted.gene])
  raw.wanted <- raw.loftee.afr[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.loftee.afr <- c(snps.with.carriers.loftee.afr, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.loftee.afr[wanted.gene] <- carriers.tmp$carrier[match(carriers.loftee.afr$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.loftee.afr)
## 18
# carriers.loftee.afr$All_Genes <- ifelse(rowSums(carriers.loftee.afr[ , -1], na.rm = TRUE) > 0, 1, 0)

#################################
## 3.a Extract snpeff European ##
#################################
raw.snpeff.eur <- raw.eur[which(colnames(raw.eur) %in% snpeff.eur$SNP)]
snpeff.eur <- snpeff.eur[snpeff.eur$SNP %in% colnames(raw.eur),]

genes <- unique(snpeff.eur$new_GENE.snpeff)
carriers.snpeff.eur <- setNames(as.data.frame(rownames(raw.eur)), "IID")
snps.with.carriers.snpeff.eur <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(snpeff.eur$SNP[snpeff.eur$new_GENE.snpeff == wanted.gene])
  raw.wanted <- raw.snpeff.eur[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.snpeff.eur <- c(snps.with.carriers.snpeff.eur, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.snpeff.eur[wanted.gene] <- carriers.tmp$carrier[match(carriers.snpeff.eur$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.snpeff.eur)
# 142
# carriers.snpeff.eur$All_Genes <- ifelse(rowSums(carriers.snpeff.eur[ , -1], na.rm = TRUE) > 0, 1, 0)

################################
## 3.b Extract snpeff African ##
################################
raw.snpeff.afr <- raw.afr[which(colnames(raw.afr) %in% snpeff.afr$SNP)]
snpeff.afr <- snpeff.afr[snpeff.afr$SNP %in% colnames(raw.afr),]

genes <- unique(snpeff.afr$new_GENE.snpeff)
carriers.snpeff.afr <- setNames(as.data.frame(rownames(raw.afr)), "IID")
snps.with.carriers.snpeff.afr <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(snpeff.afr$SNP[snpeff.afr$new_GENE.snpeff == wanted.gene])
  raw.wanted <- raw.snpeff.afr[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.snpeff.afr <- c(snps.with.carriers.snpeff.afr, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.snpeff.afr[wanted.gene] <- carriers.tmp$carrier[match(carriers.snpeff.afr$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.snpeff.afr)
## 13
# carriers.snpeff.afr$All_Genes <- ifelse(rowSums(carriers.snpeff.afr[ , -1], na.rm = TRUE) > 0, 1, 0)



# Save each object as an RDS file in the specified path
saveRDS(carriers.clinvar.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al_POI_genetics_carriers_clinvar_eur.rds")
saveRDS(carriers.clinvar.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al_POI_genetics_carriers_clinvar_afr.rds")
saveRDS(carriers.loftee.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al_POI_genetics_carriers_loftee_eur.rds")
saveRDS(carriers.loftee.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al_POI_genetics_carriers_loftee_afr.rds")
saveRDS(carriers.snpeff.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al_POI_genetics_carriers_snpeff_eur.rds")
saveRDS(carriers.snpeff.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al_POI_genetics_carriers_snpeff_afr.rds")
