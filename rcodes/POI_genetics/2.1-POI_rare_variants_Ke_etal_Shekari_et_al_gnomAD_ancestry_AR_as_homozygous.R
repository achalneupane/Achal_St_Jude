rm(list=ls())
## read WES annotation files for POI option 1
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/")

## read Ke et al genes
Ke_genes <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_et_al/Ke_et_al_genes.txt", header= T)
Shekari_genes <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Shekari_et_al/Shekari_et_al_genes.txt", header= T, sep = "\t")
Shekari_genes.AR <- Shekari_genes[Shekari_genes$Reported.gene.inheritance=="AR",] ## AR genes in Shekari

Ke_Shekari_genes <- as.data.frame(unique(c(Ke_genes$Gene, Shekari_genes$Gene)))
colnames(Ke_Shekari_genes) <- "Gene"

Ke_Shekari_genes$AR_based_on_shekari <- ifelse(Ke_Shekari_genes$Gene %in% Shekari_genes.AR$Gene, "Yes", "No")

## 1. clinvar
# clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC//annotation//snpEff_round3_preQC/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
# 155798    216
length(unique(clinvar$ANN....GENE))
# [1] 4844

clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)

tt <- cbind.data.frame(clinvar$SNP, clinvar$ID, clinvar$AF_nfe, clinvar$AF_afr, clinvar$AF)

# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar$new_GENE.clinvar <- gsub("-.*", "", clinvar$ANN....GENE)
clinvar <- cbind.data.frame(CHROM=clinvar$CHROM, POS=clinvar$POS, REF=clinvar$REF, ALT=clinvar$ALT, SNP=clinvar$SNP, ID=clinvar$ID, GENE=clinvar$ANN....GENE, new_GENE.clinvar=clinvar$new_GENE.clinvar, AF_nfe=clinvar$AF_nfe, AF_afr=clinvar$AF_afr, AF=clinvar$AF, Effect=clinvar$ANN....EFFECT, CLNSIG=clinvar$CLNSIG)
clinvar.save <- clinvar

Ke_Shekari_genes$Gene[!Ke_Shekari_genes$Gene %in% clinvar$new_GENE.clinvar]
# [1] "ANKRD31"   "ATG9A"     "BMP15"     "BMPR1A"    "BMPR1B"    "BNC1"      "C14orf39"  "DIAPH2"    "DMC1"      "EIF4ENIF1" "ESR2"      "FIGLA"     "FOXE1"    
# [14] "KHDRBS1"   "MCM9"      "NANOS3"    "NOTCH2"    "PGRMC1"    "POF1B"     "PRDM9"     "PSMC3IP"   "RAD51"     "SGO2"      "SYCP2L"    "XRCC2"     "CDKN1B"   
# [27] "DACH2"     "ERAL1"     "FMN2"      "HROB"      "IGSF10"    "LHX8"      "POU5F1"    "REC8"      "SMC1B"     "SOHLH2"    "YTHDC2" 


## Keep only those in Ke Shekari et al
clinvar <- clinvar[clinvar$new_GENE.clinvar %in% Ke_Shekari_genes$Gene,]
dim(clinvar)
# 7698   13

clinvar.all <- clinvar[which(clinvar$AF <0.01),]
clinvar.eur <- clinvar[which(clinvar$AF_nfe < 0.01),]
clinvar.afr <- clinvar[which(clinvar$AF_afr < 0.01),]




loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC//annotation/snpEff_round3_preQC//loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
# [1] 143998     90

loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)
# - 
# 139452   


# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
loftee$new_GENE.loftee <- gsub("-.*", "", loftee$SYMBOL)

Ke_Shekari_genes$Gene[!Ke_Shekari_genes$Gene %in% loftee$new_GENE.loftee]
# [1] "ATG9A"     "BMP15"     "BMPR1A"    "CYP17A1"   "EIF2B5"    "EIF4ENIF1" "FIGLA"     "FOXE1"     "FOXL2"     "GALT"      "MCM9"      "NANOS3"    "NOBOX"    
# [14] "NOG"       "NR5A1"     "NSMCE2"    "PGRMC1"    "POLR3H"    "SALL4"     "XRCC2"     "AR"        "CDKN1B"    "CPEB1"     "ERAL1"     "FMN2"      "FOXO4"    
# [27] "GGPS1"     "POU5F1"    "ZSWIM7"

## Keep only those in ACMG
loftee <- loftee[loftee$new_GENE.loftee %in% Ke_Shekari_genes$Gene,]


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

Ke_Shekari_genes$Gene[!Ke_Shekari_genes$Gene %in% snpeff$new_GENE.snpeff]
# [1] "ANKRD31"   "BLM"       "BRCA2"     "C14orf39"  "CYP17A1"   "DMC1"      "EIF4ENIF1" "FANCC"     "FANCG"     "FANCL"     "FANCM"     "FMR1"      "FOXL2"    
# [14] "FSHR"      "HAX1"      "HFM1"      "HSF2BP"    "KHDRBS1"   "LRPPRC"    "MCM8"      "MCM9"      "MEIOB"     "MSH5"      "NANOS3"    "NOBOX"     "NOG"      
# [27] "NSMCE2"    "NUP107"    "PGRMC1"    "POF1B"     "POR"       "PRDM9"     "PSMC3IP"   "RAD51"     "RECQL4"    "SALL4"     "SGO2"      "SOHLH1"    "SPATA22"  
# [40] "SPIDR"     "SYCE1"     "SYCP2L"    "WDR62"     "WT1"       "AR"        "CPEB1"     "ERAL1"     "FMN2"      "GGPS1"     "HROB"      "IGSF10"    "KASH5"    
# [53] "NHEJ1"     "POU5F1"    "REC8"      "SMC1B"     "SOHLH2"    "TRIM37"    "YTHDC2"    "ZSWIM7" 


## Keep only those in Ke and Shekari et al
snpeff <- snpeff[snpeff$new_GENE.snpeff %in% Ke_Shekari_genes$Gene,]

# make rare
snpeff.all <- snpeff[which(snpeff$AF <0.01),]
snpeff.eur <- snpeff[which(snpeff$AF_nfe < 0.01),]
snpeff.afr <- snpeff[which(snpeff$AF_afr < 0.01),]


length(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
# 1158
cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
# 1158
colnames(cc) <- "SNP"

## QCed Bim
library(data.table)
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES_QC/biallelic2//plink_all/chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated.bim")
cc$SNP[!cc$SNP %in% bim.QC$V2]
cc.final <- cc[cc$SNP %in% bim.QC$V2,]

write.table(cc.final, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_et_al_POI_rare_variants_to_extract.txt", row.names = F, col.names = F, quote = F)
length(cc.final)
# 656

## Autosomal recessive gene variants
AR.genes <- Ke_Shekari_genes$Gene[Ke_Shekari_genes$AR_based_on_shekari == "Yes"]
AR.genes.clinvar <- clinvar$SNP[clinvar$new_GENE.clinvar %in% AR.genes]
AR.genes.loftee <- loftee$SNP[loftee$new_GENE.loftee %in% AR.genes]
AR.genes.snpeff <- snpeff$SNP[snpeff$new_GENE.snpeff %in% AR.genes]
AR.variants <- (unique(c(AR.genes.clinvar, AR.genes.loftee, AR.genes.snpeff)))


## QCed for GQ, DP and VQSR only
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/plink_all/chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.bim")
table(cc$SNP %in% bim.QC$V2)
# FALSE  TRUE 
# 480   678 

## Final SJLIFE QCed data
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES_QC/biallelic2/plink_all/chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated.bim")
table(cc$SNP %in% bim.QC$V2)
# FALSE  TRUE 
# 502   656 


## Create carrier status
raw <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/POI_genetics//Ke_Shekari_et_al_POI_rare_variants_ALL_recodeA.raw")
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
table(colnames(raw) %in% cc$SNP)
# TRUE 
# 656
## Looks good!

# Make AR gene variants either 0 or 2
AR.variants.raw <- colnames(raw)[colnames(raw) %in% AR.variants]
raw[, AR.variants.raw] <- ifelse(raw[, AR.variants.raw] == 2, 2, 0)


## Admixture classification
admixture <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/admixture.rds")
EUR <- admixture$INDIVIDUAL[admixture$ancestry =="EUR"]
AFR <- admixture$INDIVIDUAL[admixture$ancestry =="AFR"]
Other <- admixture$INDIVIDUAL[admixture$ancestry == "Other"]



# table(EUR.sjlife$V2 %in% EUR$IID[EUR$cohort=="sjlife"])
# table(AFR.sjlife$V2 %in% AFR$IID[AFR$cohort=="sjlife"])

table(rownames(raw) %in% EUR)
# FALSE  TRUE
# 7761  5965 

## Extract European and African genotype data
raw.eur <- raw[rownames(raw) %in% EUR,]
raw.afr <- raw[rownames(raw) %in% AFR,]
raw.all <- raw[rownames(raw) %in% Other,]
dim(raw.eur)
# [1] 5965  656
dim(raw.afr)
# [1] 801 656
dim(raw.all)
# 551 656
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
# 69
dim(carriers.clinvar.eur)
# [1] 5965   21
## Adding carrier status for "All_Genes"
## Check if any gene column in carriers.clinvar.eur has carrier status (1) across all genes
# carriers.clinvar.eur$All_Genes <- ifelse(rowSums(carriers.clinvar.eur[ , -1], na.rm = TRUE) > 0, 1, 0)
dim(carriers.clinvar.eur)
# 5965   21
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
## 7
# carriers.clinvar.afr$All_Genes <- ifelse(rowSums(carriers.clinvar.afr[ , -1], na.rm = TRUE) > 0, 1, 0)

###############################
## 1.c Extract clinvar Other ##
###############################
raw.clinvar.all <- raw.all[which(colnames(raw.all) %in% clinvar.all$SNP)]
clinvar.all <- clinvar.all[clinvar.all$SNP %in% colnames(raw.all),]

genes <- unique(clinvar.all$new_GENE.clinvar)
carriers.clinvar.all <- setNames(as.data.frame(rownames(raw.all)), "IID")
snps.with.carriers.clinvar.all <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(clinvar.all$SNP[clinvar.all$new_GENE.clinvar == wanted.gene])
  raw.wanted <- raw.clinvar.all[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.clinvar.all <- c(snps.with.carriers.clinvar.all, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.clinvar.all[wanted.gene] <- carriers.tmp$carrier[match(carriers.clinvar.all$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.clinvar.all)
# 13
dim(carriers.clinvar.all)
# [1] 5965   10
## Adding carrier status for "All_Genes"
## Check if any gene column in carriers.clinvar.all has carrier status (1) across all genes
# carriers.clinvar.all$All_Genes <- ifelse(rowSums(carriers.clinvar.all[ , -1], na.rm = TRUE) > 0, 1, 0)
dim(carriers.clinvar.all)

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
## 20
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
## 4
# carriers.loftee.afr$All_Genes <- ifelse(rowSums(carriers.loftee.afr[ , -1], na.rm = TRUE) > 0, 1, 0)

##############################
## 2.c Extract loftee Other ##
##############################
raw.loftee.all <- raw.all[which(colnames(raw.all) %in% loftee.all$SNP)]
loftee.all <- loftee.all[loftee.all$SNP %in% colnames(raw.all),]

genes <- unique(loftee.all$new_GENE.loftee)
carriers.loftee.all <- setNames(as.data.frame(rownames(raw.all)), "IID")
snps.with.carriers.loftee.all <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(loftee.all$SNP[loftee.all$new_GENE.loftee == wanted.gene])
  raw.wanted <- raw.loftee.all[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.loftee.all <- c(snps.with.carriers.loftee.all, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.loftee.all[wanted.gene] <- carriers.tmp$carrier[match(carriers.loftee.all$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.loftee.all)
## 2
# carriers.loftee.all$All_Genes <- ifelse(rowSums(carriers.loftee.all[ , -1], na.rm = TRUE) > 0, 1, 0)


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
# 82
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
## 7
# carriers.snpeff.afr$All_Genes <- ifelse(rowSums(carriers.snpeff.afr[ , -1], na.rm = TRUE) > 0, 1, 0)


##############################
## 3.c Extract snpeff Other ##
##############################
raw.snpeff.all <- raw.all[which(colnames(raw.all) %in% snpeff.all$SNP)]
snpeff.all <- snpeff.all[snpeff.all$SNP %in% colnames(raw.all),]

genes <- unique(snpeff.all$new_GENE.snpeff)
carriers.snpeff.all <- setNames(as.data.frame(rownames(raw.all)), "IID")
snps.with.carriers.snpeff.all <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(snpeff.all$SNP[snpeff.all$new_GENE.snpeff == wanted.gene])
  raw.wanted <- raw.snpeff.all[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.snpeff.all <- c(snps.with.carriers.snpeff.all, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.snpeff.all[wanted.gene] <- carriers.tmp$carrier[match(carriers.snpeff.all$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.snpeff.all)
# 12
# carriers.snpeff.all$All_Genes <- ifelse(rowSums(carriers.snpeff.all[ , -1], na.rm = TRUE) > 0, 1, 0)



# Save each object as an RDS file in the specified path
saveRDS(carriers.clinvar.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_clinvar_gnomAD_eur_AR_as_homozygous.rds")
saveRDS(carriers.clinvar.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_clinvar_gnomAD_afr_AR_as_homozygous.rds")
saveRDS(carriers.clinvar.all, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_clinvar_gnomAD_other_AR_as_homozygous.rds")

saveRDS(carriers.loftee.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_loftee_gnomAD_eur_AR_as_homozygous.rds")
saveRDS(carriers.loftee.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_loftee_gnomAD_afr_AR_as_homozygous.rds")
saveRDS(carriers.loftee.all, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_loftee_gnomAD_other_AR_as_homozygous.rds")

saveRDS(carriers.snpeff.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_snpeff_gnomAD_eur_AR_as_homozygous.rds")
saveRDS(carriers.snpeff.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_snpeff_gnomAD_afr_AR_as_homozygous.rds")
saveRDS(carriers.snpeff.all, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ke_Shekari_etal_POI_genetics_carriers_snpeff_gnomAD_other_AR_as_homozygous.rds")



