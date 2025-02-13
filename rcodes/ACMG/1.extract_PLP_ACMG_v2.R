## This is the updated analysis based on Cindy's request:
## See Lauren's email on 2/5/2025 and Cindy's email on 2/3/2025
# The possible categories are Cancer, Cardiovascular, Cardiovascular Metabolic, Metabolic and Miscellaneous. 
# I combined Cardiovascular with Cardiovascular Metabolic into Cardiovascular. I combined Metabolic with Miscellaneous into Miscellaneous. 

# Disease category cardiovascular included genes with cardiovascular and
# cardiovascular metabolic phenotypes. Disease category miscellaneous included
# genes with miscellaneous and metabolic phenotypes. Disease phenotype cancer
# only included genes with cancer phenotypes.
library(data.table)
rm(list=ls())
## Read ACMG genes
# ACMG <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/ACMG_genes.txt", header = T, sep = "\t", stringsAsFactors = F)
ACMG <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/ACMG2.txt", header = T, sep = "\t", stringsAsFactors = F) ## This is from Lauren
# ACMG$group <- ACMG$Disease_name_and_MIM_number
# ACMG$group <- sub("\\(.*", "", ACMG$group)
# ACMG$group <- sub(",.*", "", ACMG$group)
# ACMG$group <- sub("[0-9].*", "", ACMG$group)
# ACMG$group <- trimws(ACMG$group)
# ACMG$group[ACMG$group == "RPE"] <- "RPE65-related retinopathy"
# ACMG$group[ACMG$group == "Juvenile polyposis/hereditary hemorrhagic telangiectasia syndrome"] <- "Juvenile polyposis syndrome"
# ACMG$group <- sub(" type*", "", ACMG$group)
# as.data.frame(table(ACMG$group))
# ACMG$GENE <- sub("<.*", "", ACMG$Gene_via_GTR)
# ACMG$inheritence <- ""
# 
# ACMG$disease <- gsub("\'", "", ACMG$group)
# ACMG$disease <- gsub(" ",".", ACMG$disease)
# disease <- unique(ACMG$disease)

# overall <- ACMG$GENE

cancer <- unique(ACMG$GENE[grepl("Cancer", ACMG$PhenotypeCategory, ignore.case = T)])
cardiovascular <- unique(ACMG$GENE[grepl("Cardiovascular", ACMG$PhenotypeCategory, ignore.case = T)])
# metabolic <- unique(ACMG$GENE[grepl("^Metabo", ACMG$PhenotypeCategory, ignore.case = T)])
# cardiovascular_and_metabolic <- unique(c(cardiovascular, metabolic))
miscellaneous <- unique(ACMG$GENE[grepl("^Miscell|^Metabolic", ACMG$PhenotypeCategory, ignore.case = T)])
# miscellaneous_and_metabolic <- unique(c(miscellaneous, metabolic))

# disease_groups <- list(cancer=cancer, cardiovascular=cardiovascular, cardiovascular_and_metabolic=cardiovascular_and_metabolic, miscellaneous_and_metabolic=miscellaneous_and_metabolic, miscellaneous=miscellaneous, metabolic=metabolic)
disease_groups <- list(cancer=cancer, cardiovascular=cardiovascular, miscellaneous=miscellaneous)



## AD or AR
# ACMG$inheritence[grepl("Adenomatous polyposis coli|Aortic aneurysm|Arrhythmogenic right ventricular cardiomyopathy|Breast-ovarian cancer|Brugada syndrome|Catecholaminergic|Dilated cardiomyopathy|Ehlers|Fabry|Familial hypercholesterolemia|Familial hypertrophic cardiomyopathy|Familial medullary thyroid carcinoma|Hereditary breast cancer|hemochromatosis|Hereditary hemorrhagic telangiectasia|Hereditary paraganglioma-pheochromocytoma|Hereditary transthyretin-related amyloidosis|Hypercholesterolemia|Juvenile polyposis|Li-Fraumeni syndrome|Loeys-Dietz syndrome|Long QT|Lynch|Malignant hyperthermia|Marfan|Maturity-Onset of Diabetes|Multiple endocrine neoplasia|Myofibrillar myopathy|Neurofibromatosis|Ornithine carbamoyltransferase deficiency|Paragangliomas|Peutz-Jeghers syndrome|Pheochromocytoma|PTEN hamartoma|Retinoblastoma|Tuberous sclerosis|Von Hippel-Lindau|Wilms", ACMG$group, ignore.case = T)] <- "AD"
# ACMG$inheritence[grepl("Biotinidase deficiency|MYH-associated polyposis|Pompe disease|RPE65|Wilson disease", ACMG$group, ignore.case = T)] <- "AR"
ACMG$inheritence <- ACMG$inheritance

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/")

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

ACMG$GENE[!ACMG$GENE %in% clinvar$new_GENE.clinvar]
## Keep only those in ACMG
clinvar <- clinvar[clinvar$new_GENE.clinvar %in% ACMG$GENE,]

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

ACMG$GENE[!ACMG$GENE %in% loftee$new_GENE.loftee]
## Keep only those in ACMG
loftee <- loftee[loftee$new_GENE.loftee %in% ACMG$GENE,]


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

ACMG$GENE[!ACMG$GENE %in% snpeff$new_GENE.snpeff]
## Keep only those in ACMG
snpeff <- snpeff[snpeff$new_GENE.snpeff %in% ACMG$GENE,]

# make rare
snpeff.all <- snpeff[which(snpeff$AF <0.01),]
snpeff.eur <- snpeff[which(snpeff$AF_nfe < 0.01),]
snpeff.afr <- snpeff[which(snpeff$AF_afr < 0.01),]



union.eur <- rbind.data.frame(cbind.data.frame(SNP=clinvar.eur$SNP, new_GENE.union=clinvar.eur$new_GENE.clinvar),
                                      cbind.data.frame(SNP=loftee.eur$SNP, new_GENE.union=loftee.eur$new_GENE.loftee),
                                      cbind.data.frame(SNP=snpeff.eur$SNP, new_GENE.union=snpeff.eur$new_GENE.snpeff))

union.eur <- union.eur[!duplicated(union.eur$SNP),]


union.afr <- rbind.data.frame(cbind.data.frame(SNP=clinvar.afr$SNP, new_GENE.union=clinvar.afr$new_GENE.clinvar),
                                      cbind.data.frame(SNP=loftee.afr$SNP, new_GENE.union=loftee.afr$new_GENE.loftee),
                                      cbind.data.frame(SNP=snpeff.afr$SNP, new_GENE.union=snpeff.afr$new_GENE.snpeff))

union.afr <- union.afr[!duplicated(union.afr$SNP),]


union.all <- rbind.data.frame(cbind.data.frame(SNP=clinvar.all$SNP, new_GENE.union=clinvar.all$new_GENE.clinvar),
                              cbind.data.frame(SNP=loftee.all$SNP, new_GENE.union=loftee.all$new_GENE.loftee),
                              cbind.data.frame(SNP=snpeff.all$SNP, new_GENE.union=snpeff.all$new_GENE.snpeff))

union.all <- union.all[!duplicated(union.all$SNP),]


length(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
# 930
cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
# 930
colnames(cc) <- "SNP"

## QCed Bim
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated.bim")
cc$SNP[!cc$SNP %in% bim.QC$V2]
cc.final <- cc[cc$SNP %in% bim.QC$V2,]

# write.table(cc.final, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/ACMG_rare_variants_to_extract.txt", row.names = F, col.names = F, quote = F)
length(cc.final)
# 728

## Autosomal recessive gene variants
AR.genes <- ACMG$GENE[ACMG$inheritence == "AR"]
AR.genes.clinvar <- clinvar$SNP[clinvar$new_GENE.clinvar %in% AR.genes]
AR.genes.loftee <- loftee$SNP[loftee$new_GENE.loftee %in% AR.genes]
AR.genes.snpeff <- snpeff$SNP[snpeff$new_GENE.snpeff %in% AR.genes]
AR.variants <- (unique(c(AR.genes.clinvar, AR.genes.loftee, AR.genes.snpeff)))

## preQC
bim.preQC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim")
cc <- as.character(cc$SNP)
cc[!cc %in% bim.preQC$V2]
table(cc %in% bim.preQC$V2)
# TRUE 
# 930

## QCed for GQ, DP and VQSR
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim")
table(cc %in% bim.QC$V2)
# FALSE  TRUE 
# 179   751

## Final SJLIFE QCed data
bim.QC.WES <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated.bim")
table(cc %in% bim.QC.WES$V2)
# FALSE  TRUE 
# 202   728 


## Create carrier status
raw <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/ACMG/ACMG_rare_variants_ALL_recodeA.raw")
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
# 728
## Looks good!

# Make AR gene variants either 0 or 2
AR.variants.raw <- colnames(raw)[colnames(raw) %in% AR.variants]
raw[, AR.variants.raw] <- ifelse(raw[, AR.variants.raw] == 2, 2, 0)


# save.image("rare_varaints_ACMG_data.RData")
# load("rare_varaints_ACMG_data.RData")

# Some ccss/sjlife overlapping samples in WGS data were renames as SJLIFE, but not in WES, so naming them back to CCSS.
sjlife_ccss_exp_overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sjlife_ccss_exp_overlaps.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/pca/Survivor_WGS_EUR_based_on_1kGP_Phase_3_data.txt", header = T)
# EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/pca/Survivor_WGS_and_1kGP_final_EUR_top_20_PCs.eigenvec", header = F)
# colnames(EUR) <- c("FID", "IID", paste0("PC",1:20))
# EUR <- EUR[!EUR$FID %in% "SJL2532507",]
## Find matching indices for SJLID in EUR$IID and replace with CCSSID
matches <- match(EUR$IID, sjlife_ccss_exp_overlaps$SJLID)
EUR$IID[!is.na(matches)] <- sjlife_ccss_exp_overlaps$CCSSID[matches[!is.na(matches)]]

AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/pca/Survivor_WGS_AFR_based_on_1kGP_Phase_3_data.txt", header = T)
# AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/pca/Survivor_WGS_and_1kGP_final_AFR_top_20_PCs.eigenvec", header = F)
# colnames(AFR) <- c("FID", "IID", paste0("PC",1:20))
## Find matching indices for SJLID in EUR$IID and replace with CCSSID
matches <- match(AFR$IID, sjlife_ccss_exp_overlaps$SJLID)
AFR$IID[!is.na(matches)] <- sjlife_ccss_exp_overlaps$CCSSID[matches[!is.na(matches)]]

sjlife_ccss_exp_overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sjlife_ccss_exp_overlaps.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ALL <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/pca/Survivor_WGS_and_1kGP_final_all_Survivor_WGS_samples_20_PCs.eigenvec", header = F)
# EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/pca/Survivor_WGS_and_1kGP_final_EUR_top_20_PCs.eigenvec", header = F)
colnames(ALL) <- c("FID", "IID", paste0("PC",1:20))
ALL <- ALL[!ALL$FID %in% "SJL2532507",]
## Find matching indices for SJLID in EUR$IID and replace with CCSSID
matches <- match(ALL$IID, sjlife_ccss_exp_overlaps$SJLID)
ALL$IID[!is.na(matches)] <- sjlife_ccss_exp_overlaps$CCSSID[matches[!is.na(matches)]]


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

# save.image("rare_variant_ACMG_group_data.RData")

# load("rare_variant_ACMG_group_data.RData")
## Extract European and African genotype data
raw.eur <- raw[rownames(raw) %in% EUR$IID,]
raw.afr <- raw[rownames(raw) %in% AFR$IID,]
afr.eur <- rownames(rbind.data.frame(raw.eur, raw.afr))
raw.all <- raw[!rownames(raw) %in% afr.eur,] # other

dim(raw.eur)
# [1] 6001  728
dim(raw.afr)
# [1] 846 728
dim(raw.all)
# 1173  728
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
# 162
dim(carriers.clinvar.eur)
# [1] 6001   44
# Adding carrier status for "All_Genes"
# Check if any gene column in carriers.clinvar.eur has carrier status (1) across all genes
carriers.clinvar.eur$All_Genes <- ifelse(rowSums(carriers.clinvar.eur[ , -1], na.rm = TRUE) > 0, 1, 0)
dim(carriers.clinvar.eur)
# [1] 6001   45

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.clinvar.eur) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.clinvar.eur[colnames(carriers.clinvar.eur) %in% wanted.genes]
#   carriers.clinvar.eur[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.clinvar.eur) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.clinvar.eur[colnames(carriers.clinvar.eur) %in% wanted.genes]
  carriers.clinvar.eur[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.clinvar.eur$cohort <- filtered_wes_samples$Source[match(carriers.clinvar.eur$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.clinvar.eur[-1])


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
# 24
carriers.clinvar.afr$All_Genes <- ifelse(rowSums(carriers.clinvar.afr[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.clinvar.afr) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.clinvar.afr[colnames(carriers.clinvar.afr) %in% wanted.genes]
#   carriers.clinvar.afr[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.clinvar.afr) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.clinvar.afr[colnames(carriers.clinvar.afr) %in% wanted.genes]
  carriers.clinvar.afr[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.clinvar.afr$cohort <- filtered_wes_samples$Source[match(carriers.clinvar.afr$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.clinvar.afr[-1])


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
# 24
carriers.clinvar.all$All_Genes <- ifelse(rowSums(carriers.clinvar.all[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.clinvar.all) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.clinvar.all[colnames(carriers.clinvar.all) %in% wanted.genes]
#   carriers.clinvar.all[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.clinvar.all) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.clinvar.all[colnames(carriers.clinvar.all) %in% wanted.genes]
  carriers.clinvar.all[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.clinvar.all$cohort <- filtered_wes_samples$Source[match(carriers.clinvar.all$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.clinvar.all[-1])


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
# 35
carriers.loftee.eur$All_Genes <- ifelse(rowSums(carriers.loftee.eur[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.loftee.eur) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.loftee.eur[colnames(carriers.loftee.eur) %in% wanted.genes]
#   carriers.loftee.eur[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.loftee.eur) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.loftee.eur[colnames(carriers.loftee.eur) %in% wanted.genes]
  carriers.loftee.eur[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.loftee.eur$cohort <- filtered_wes_samples$Source[match(carriers.loftee.eur$IID, filtered_wes_samples$SampleID)]

result_table <- count_zeros_ones(carriers.loftee.eur[-1])

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
# 8
carriers.loftee.afr$All_Genes <- ifelse(rowSums(carriers.loftee.afr[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.loftee.afr) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.loftee.afr[colnames(carriers.loftee.afr) %in% wanted.genes]
#   carriers.loftee.afr[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.loftee.afr) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.loftee.afr[colnames(carriers.loftee.afr) %in% wanted.genes]
  carriers.loftee.afr[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.loftee.afr$cohort <- filtered_wes_samples$Source[match(carriers.loftee.afr$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.loftee.afr[-1])


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
# 8
carriers.loftee.all$All_Genes <- ifelse(rowSums(carriers.loftee.all[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.loftee.all) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.loftee.all[colnames(carriers.loftee.all) %in% wanted.genes]
#   carriers.loftee.all[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.loftee.all) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.loftee.all[colnames(carriers.loftee.all) %in% wanted.genes]
  carriers.loftee.all[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.loftee.all$cohort <- filtered_wes_samples$Source[match(carriers.loftee.all$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.loftee.all[-1])


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
# 203
carriers.snpeff.eur$All_Genes <- ifelse(rowSums(carriers.snpeff.eur[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.snpeff.eur) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.snpeff.eur[colnames(carriers.snpeff.eur) %in% wanted.genes]
#   carriers.snpeff.eur[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.snpeff.eur) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.snpeff.eur[colnames(carriers.snpeff.eur) %in% wanted.genes]
  carriers.snpeff.eur[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.snpeff.eur$cohort <- filtered_wes_samples$Source[match(carriers.snpeff.eur$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.snpeff.eur[-1])

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
# 33
carriers.snpeff.afr$All_Genes <- ifelse(rowSums(carriers.snpeff.afr[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.snpeff.afr) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.snpeff.afr[colnames(carriers.snpeff.afr) %in% wanted.genes]
#   carriers.snpeff.afr[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.snpeff.afr) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.snpeff.afr[colnames(carriers.snpeff.afr) %in% wanted.genes]
  carriers.snpeff.afr[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.snpeff.afr$cohort <- filtered_wes_samples$Source[match(carriers.snpeff.afr$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.snpeff.afr[-1])

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
# 33
carriers.snpeff.all$All_Genes <- ifelse(rowSums(carriers.snpeff.all[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.snpeff.all) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.snpeff.all[colnames(carriers.snpeff.all) %in% wanted.genes]
#   carriers.snpeff.all[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.snpeff.all) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.snpeff.all[colnames(carriers.snpeff.all) %in% wanted.genes]
  carriers.snpeff.all[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.snpeff.all$cohort <- filtered_wes_samples$Source[match(carriers.snpeff.all$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.snpeff.all[-1])


###############################
## Clinvar, Loftee or SnpEff ##
###############################
####################################################
## 4.a Extract Clinvar, Loftee or SnpEff European ##
####################################################
raw.union.eur <- raw.eur[which(colnames(raw.eur) %in% union.eur$SNP)]
union.eur <- union.eur[union.eur$SNP %in% colnames(raw.eur),]

genes <- unique(union.eur$new_GENE.union)
carriers.union.eur <- setNames(as.data.frame(rownames(raw.eur)), "IID")
snps.with.carriers.union.eur <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(union.eur$SNP[union.eur$new_GENE.union == wanted.gene])
  raw.wanted <- raw.union.eur[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.union.eur <- c(snps.with.carriers.union.eur, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.union.eur[wanted.gene] <- carriers.tmp$carrier[match(carriers.union.eur$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.union.eur)
# 384
carriers.union.eur$All_Genes <- ifelse(rowSums(carriers.union.eur[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.union.eur) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.union.eur[colnames(carriers.union.eur) %in% wanted.genes]
#   carriers.union.eur[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.union.eur) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.union.eur[colnames(carriers.union.eur) %in% wanted.genes]
  carriers.union.eur[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.union.eur$cohort <- filtered_wes_samples$Source[match(carriers.union.eur$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.union.eur[-1])

###################################################
## 4.b Extract Clinvar, Loftee or SnpEff African ##
###################################################
raw.union.afr <- raw.afr[which(colnames(raw.afr) %in% union.afr$SNP)]
union.afr <- union.afr[union.afr$SNP %in% colnames(raw.afr),]

genes <- unique(union.afr$new_GENE.union)
carriers.union.afr <- setNames(as.data.frame(rownames(raw.afr)), "IID")
snps.with.carriers.union.afr <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(union.afr$SNP[union.afr$new_GENE.union == wanted.gene])
  raw.wanted <- raw.union.afr[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.union.afr <- c(snps.with.carriers.union.afr, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.union.afr[wanted.gene] <- carriers.tmp$carrier[match(carriers.union.afr$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.union.afr)
# 33
carriers.union.afr$All_Genes <- ifelse(rowSums(carriers.union.afr[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.union.afr) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.union.afr[colnames(carriers.union.afr) %in% wanted.genes]
#   carriers.union.afr[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.union.afr) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.union.afr[colnames(carriers.union.afr) %in% wanted.genes]
  carriers.union.afr[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.union.afr$cohort <- filtered_wes_samples$Source[match(carriers.union.afr$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.union.afr[-1])


#################################################
## 4.b Extract Clinvar, Loftee or SnpEff Other ##
#################################################
raw.union.all <- raw.all[which(colnames(raw.all) %in% union.all$SNP)]
union.all <- union.all[union.all$SNP %in% colnames(raw.all),]

genes <- unique(union.all$new_GENE.union)
carriers.union.all <- setNames(as.data.frame(rownames(raw.all)), "IID")
snps.with.carriers.union.all <- c()

## Status per gene
for (i in 1:length(genes)){
  wanted.gene <- genes[i]
  print(paste0("Doing gene ", wanted.gene))
  wanted.snps <- unique(union.all$SNP[union.all$new_GENE.union == wanted.gene])
  raw.wanted <- raw.union.all[wanted.snps]
  snps.with.carriers.tmp <- names(colSums(raw.wanted, na.rm = T))[colSums(raw.wanted, na.rm = T) > 0]
  if(length(snps.with.carriers.tmp) == 0){
    next
  }
  raw.wanted <- raw.wanted[snps.with.carriers.tmp]
  snps.with.carriers.union.all <- c(snps.with.carriers.union.all, snps.with.carriers.tmp)
  carriers.tmp <- setNames(as.data.frame(rownames(raw.wanted)), "IID")
  carriers.tmp$carrier <- ifelse(rowSums(raw.wanted[grepl("chr", colnames(raw.wanted))], na.rm = T) > 0, 1, 0)
  carriers.union.all[wanted.gene] <- carriers.tmp$carrier[match(carriers.union.all$IID, carriers.tmp$IID)]
}

length(snps.with.carriers.union.all)
# 33
carriers.union.all$All_Genes <- ifelse(rowSums(carriers.union.all[ , -1], na.rm = TRUE) > 0, 1, 0)

# ## Status per ACMG disease groups
# for (i in 1:length(disease)){
#   wanted.disease <- disease[i]
#   print(paste0("Doing disease ", wanted.disease))
#   wanted.genes <- unique(ACMG$GENE[ACMG$disease == wanted.disease])
#   if(sum(colnames(carriers.union.all) %in% wanted.genes) == 0){
#     next
#   }
#   raw.wanted <- carriers.union.all[colnames(carriers.union.all) %in% wanted.genes]
#   carriers.union.all[wanted.disease] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
# }

## disease groups per Jensson et al
for (i in 1:length(disease_groups)){
  wanted.disease <- unlist(unname(disease_groups[i]))
  print(paste0("Doing disease ", names(disease_groups[i])))
  wanted.genes <- unique(ACMG$GENE[ACMG$GENE %in% wanted.disease])
  if(sum(colnames(carriers.union.all) %in% wanted.genes) == 0){
    next
  }
  raw.wanted <- carriers.union.all[colnames(carriers.union.all) %in% wanted.genes]
  carriers.union.all[names(disease_groups[i])] <- ifelse(rowSums(raw.wanted, na.rm = T) > 0, 1, 0)
}

carriers.union.all$cohort <- filtered_wes_samples$Source[match(carriers.union.all$IID, filtered_wes_samples$SampleID)]
result_table <- count_zeros_ones(carriers.union.all[-1])



# Save each object as an RDS file in the specified path
saveRDS(carriers.clinvar.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_clinvar_eur.rds")
saveRDS(carriers.clinvar.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_clinvar_afr.rds")
saveRDS(carriers.clinvar.all, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_clinvar_other.rds")

saveRDS(carriers.loftee.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_loftee_eur.rds")
saveRDS(carriers.loftee.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_loftee_afr.rds")
saveRDS(carriers.loftee.all, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_loftee_other.rds")

saveRDS(carriers.snpeff.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_snpeff_eur.rds")
saveRDS(carriers.snpeff.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_snpeff_afr.rds")
saveRDS(carriers.snpeff.all, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_snpeff_other.rds")

saveRDS(carriers.union.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_union_eur.rds")
saveRDS(carriers.union.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_union_afr.rds")
saveRDS(carriers.union.all, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/carriers_union_other.rds")


## CSG60
# See 1. extract_PLP_CSG60_v2.R code

EUR <- EUR [!EUR$FID %in% "SJL2532507",]

table(EUR$IID %in% carriers.clinvar.eur$IID)
table(AFR$IID %in% carriers.clinvar.afr$IID)

EUR <- EUR[EUR$IID %in% carriers.clinvar.eur$IID,]
AFR <- AFR[AFR$IID %in% carriers.clinvar.afr$IID,]

saveRDS(EUR, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/EUR.PCA.rds")
saveRDS(AFR, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/AFR.PCA.rds")



ACMG <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/ACMG2.txt", header = T, sep = "\t", stringsAsFactors = F) ## From Lauren
#################################################
## Cardiovascular and Cardiovascular Metabolic ##
#################################################
cardiovascular.pheno <- ACMG$PhenotypeCategory[grepl("Cardiovascular", ACMG$PhenotypeCategory, ignore.case = T)]
table(cardiovascular.pheno) # phenotype category
# Cardiovascular Cardiovascular Metabolic 
# 52                        2 
cardiovascular <- unique(ACMG$GENE[grepl("Cardiovascular", ACMG$PhenotypeCategory, ignore.case = T)])
cardiovascular # 41 genes
# cardiovascular
# [1] "ACTA2"  "ACTC1"  "APOB"   "BAG3"   "CALM1"  "CALM2"  "CALM3"  "CASQ2"  "COL3A1" "DES"    "DSC2"   "DSG2"   "DSP"    "FBN1"  
# [15] "FLNC"   "GLA"    "KCNH2"  "KCNQ1"  "LDLR"   "LMNA"   "MYBPC3" "MYH11"  "MYH7"   "MYL2"   "MYL3"   "PCSK9"  "PKP2"   "PRKAG2"
# [29] "RBM20"  "RYR2"   "SCN5A"  "SMAD3"  "TGFBR1" "TGFBR2" "TMEM43" "TNNC1"  "TNNI3"  "TNNT2"  "TPM1"   "TRDN"   "TTN" 

############
## Cancer ##
############
cancer.pheno <- ACMG$PhenotypeCategory[grepl("cancer", ACMG$PhenotypeCategory, ignore.case = T)]
table(cancer.pheno) # phenotype category
# Cancer 
# 30
cancer <- unique(ACMG$GENE[grepl("Cancer", ACMG$PhenotypeCategory, ignore.case = T)])
cancer # 28 genes
# [1] "APC"     "BMPR1A"  "BRCA1"   "BRCA2"   "MAX"     "MEN1"    "MLH1"    "MSH2"    "MSH6"    "MUTYH"   "NF2"     "PALB2"  
# [13] "PMS2"    "PTEN"    "RB1"     "RET"     "SDHAF2"  "SDHB"    "SDHC"    "SDHD"    "SMAD4"   "STK11"   "TMEM127" "TP53"   
# [25] "TSC1"    "TSC2"    "VHL"     "WT1"

###################
## Miscellaneous ##
###################
miscellaneous.pheno <- ACMG$PhenotypeCategory[grepl("^Miscell|^Metabolic", ACMG$PhenotypeCategory, ignore.case = T)]
table(miscellaneous.pheno) # phenotype category
# Metabolic Miscellaneous 
# 3            10 

miscellaneous <- unique(ACMG$GENE[grepl("^Miscell|^Metabolic", ACMG$PhenotypeCategory, ignore.case = T)])
miscellaneous # 13 genes
# [1] "ACVRL1"  "ATP7B"   "BTD"     "CACNA1S" "ENG"     "GAA"     "HFE"     "HNF1A"   "OTC"     "RPE65"   "RYR1"    "SMAD4"  
# [13] "TTR" 

lauren <- read.delim("Mills_ACMG_collapsedPhenotypes.tsv", header = T, sep = "\t")
table(lauren$collapsedPhenotype)
lauren.cancer <- lauren$Gene[lauren$collapsedPhenotype=="Cancer"]
lauren.cancer

lauren.cardio.samples <- read.delim("Mills_cardio_ccss.tsv", header = T, sep = "\t")

####################
## carrier counts ##
####################
# carriers.clinvar.eur <- readRDS("C:/Users/aneupane/Downloads/carriers_clinvar_eur.rds")
cc <- carriers.clinvar.eur[carriers.clinvar.eur$cohort == "ccss_exp",]
cc <- carriers.union.eur[carriers.union.eur$cohort == "ccss_exp",]
dim(cc)
# 2296
table(cc$cancer) ## Seems Ok!
# 0    1 
# 2266   31
table(cc$cardiovascular) # Table has 26
# 0    1 
# 2273   25 
table(cc$miscellaneous) # Table has 15
# 0    1 
# 2282   15 

table(rowSums(cc[colnames(cc) %in% cancer]))
# 0    1 
# 2265   31
table(rowSums(cc[colnames(cc) %in% cardiovascular]))
# 0    1 
# 2271   25 
table(rowSums(cc[colnames(cc) %in% miscellaneous]))
# 0    1 
# 2281   15


gg <- cc$IID[cc$cardiovascular == 1]
lauren.cardio.samples$IID[!lauren.cardio.samples$IID %in% gg]
# gg <- readRDS("C:/Users/aneupane/Downloads/carriers_clinvar_eur.rds")
# cc <- gg[gg$cohort == "ccss_exp",]

# 24447368

######################
## Identify samples ##
######################
vcf.samples.wgs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WGS_QCed/Phenotypes//final_WGS_samples.txt", header = T)
table(ALL$IID %in% vcf.samples.wgs$VCFrename )
# FALSE  TRUE 
# 5  7829
table(vcf.samples.wgs$VCFrename %in% rownames(raw))
# FALSE  TRUE 
# 37  7942
vcf.samples.wgs$YESNO <- vcf.samples.wgs$VCFrename %in% admixture$INDIVIDUAL


## Admixture
admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
admixture$cohort [grepl("_", admixture$INDIVIDUAL)] <- "CCSS_original"
admixture$INDIVIDUAL <- sub("_.*", "", admixture$INDIVIDUAL)
admixture$ancestry <- NA
admixture$ancestry[admixture$EUR > 0.8] <- "EUR"
admixture$ancestry[admixture$AFR > 0.6] <- "AFR"
admixture$ancestry[is.na(admixture$ancestry)] <- "Other"

# Find the matching indices
matched_idx <- match(admixture$INDIVIDUAL, vcf.samples.wgs$VCFrename)
# Update only non-NA matches without overwriting existing values
update_mask <- !is.na(matched_idx) & is.na(admixture$cohort)  # Only update if cohort is NA
admixture$cohort[update_mask] <- vcf.samples.wgs$pop[matched_idx[update_mask]]

saveRDS(admixture, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/admixture.rds")
## Check
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
cc <- admixture[admixture$cohort == "Survivor",]
table(cc$INDIVIDUAL %in% PHENO.ANY_SN$sjlid)
# FALSE  TRUE 
# 2  4393 
cc[!cc$INDIVIDUAL %in% PHENO.ANY_SN$sjlid,]

sequencing.Record.ccss_exp <- read_excel("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/SJLIFE_WGS_and_CCSS_SNP_WES_WGS_samples_and_overlap.xlsx", sheet = "CCSS_exp_WGS_2839", col_names = F)
table(admixture$INDIVIDUAL %in% sequencing.Record.ccss_exp$...1)
# FALSE  TRUE 
# 10240  2839 