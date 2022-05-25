######################
## Achal Neupane    ##
## Date: 05/12/2022 ##
######################
## First check how many variants in Zhaoming et al and Qin et al are also in our VCF dataset
# Read variant lists from Zhaoming et al 
library(data.table)
library(dplyr)
Sys.setlocale("LC_ALL", "C")
zhaoming.etal.vars <- read.delim("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Zhaoming_Wang_Genetic_risk_for_subsequent_neoplasms_JCO.2018.77.8589/P-PL-sjlife-genetics-sn_SNV_INDELS.txt", sep = "\t", header =T, stringsAsFactors = F)
head(zhaoming.etal.vars)
dim(zhaoming.etal.vars)
var.classification <- table(zhaoming.etal.vars$Classification)
var.classification
# LP   P 
# 160 188

zhaoming.etal.vars$STUDY <- "Genetic_Risk_for_SN"
zhaoming.etal.vars$KEY.pos <- paste0("chr", zhaoming.etal.vars$Chr, ":", zhaoming.etal.vars$Pos_GRCh38)
zhaoming.etal.vars$KEY.varID <- paste0("chr", zhaoming.etal.vars$Chr, ":", zhaoming.etal.vars$Pos_GRCh38, ":", zhaoming.etal.vars$Reference_Allele, ":", zhaoming.etal.vars$Mutant_Allele)

## Reading variants from Qin et al (DNA-repair)
qin.etal.vars <- read.table("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Na_Qin_Pathogenic Germline Mutations in DNA Repair Genes in Combination With Cancer Treatment Exposures/supplementary_DS_jco.19.02760.txt", sep = "\t", header =T, stringsAsFactors = F)
head(qin.etal.vars)

qin.etal.vars$STUDY <- "DNA_Repair"
qin.etal.vars$KEY.pos <- paste0("chr", qin.etal.vars$Chr, ":", qin.etal.vars$Pos_GRCh38)
qin.etal.vars$KEY.varID <- paste0("chr", qin.etal.vars$Chr, ":", qin.etal.vars$Pos_GRCh38, ":", qin.etal.vars$Reference_Allele, ":", qin.etal.vars$Mutant_Allele)

## Check if the variants in Zhaoming et al and Qin et al are also in our datasets
zhaoming.etal.vars$Chr <- trimws(zhaoming.etal.vars$Chr, which = "both")
zhaoming.etal.vars$Pos_GRCh38 <- trimws(zhaoming.etal.vars$Pos_GRCh38, which = "both")

qin.etal.vars$Chr <- trimws(qin.etal.vars$Chr, which = "both")
qin.etal.vars$Pos_GRCh38 <- trimws(qin.etal.vars$Pos_GRCh38, which = "both")

## Check if these variants are in our VCF files
i=5

Qin_in_VCF <- {}
Zhaoming_in_VCF <- {}
# CHR=1:22
CHR=20
for(i in 1:length(CHR)){
print(paste0("Doing Chr ", CHR[i]))
VCF_vars <- fread(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/variants_in_VCF/", "CHR", CHR[i],"_Vars.vcf"), sep = "\t")
VCF_vars$SNP_pos <- substr(VCF_vars$ID, 1, sapply(gregexpr(":", VCF_vars$ID), "[", 2) - 1)
print(paste0("Qin ", sum(VCF_vars$SNP_pos %in% qin.etal.vars$KEY.pos)))
print(paste0("Zhaoming ", sum(VCF_vars$SNP_pos %in% zhaoming.etal.vars$KEY.pos)))
Qin_in_VCF.tmp <- VCF_vars[VCF_vars$SNP_pos %in% qin.etal.vars$KEY.pos,]
Zhaoming_in_VCF.tmp <- VCF_vars[VCF_vars$SNP_pos %in% zhaoming.etal.vars$KEY.pos,]
Qin_in_VCF <- rbind.data.frame(Qin_in_VCF, Qin_in_VCF.tmp)
Zhaoming_in_VCF <- rbind.data.frame(Zhaoming_in_VCF, Zhaoming_in_VCF.tmp)
}

write.table(Qin_in_VCF, "Qin_in_VCF.txt", sep = "\t", quote = FALSE, row.names = F, col.names =T, append = TRUE)
write.table(Zhaoming_in_VCF, "Zhaoming_in_VCF.txt", sep = "\t", quote = FALSE, row.names = F, col.names =T, append = TRUE)

length(unique(qin.etal.vars$KEY.pos))
# 389
sum(unique(Qin_in_VCF$SNP_pos) %in% unique(qin.etal.vars$KEY.pos))
# 205

length(unique(zhaoming.etal.vars$KEY.pos))
# 295
sum(unique(Zhaoming_in_VCF$ID) %in% unique(zhaoming.etal.vars$KEY.varID))
# 159


sum(Qin_in_VCF$ID %in% qin.etal.vars$KEY.varID)
# 205
Qin.types <- qin.etal.vars[qin.etal.vars$KEY.varID %in% Qin_in_VCF$ID,]
Qin_in_VCF$ID %in% qin.etal.vars$KEY.varID

sum(Zhaoming_in_VCF$ID %in% zhaoming.etal.vars$KEY.varID)
# 159



## fix positions for Indels
sum(grepl("^-$", qin.etal.vars$START))


## Creating bed files from variants in Qin and Zhaoming et al.
zhaoming.etal.vars$KEY.pos

cbind(sapply(strsplit(zhaoming.etal.vars$KEY.pos,":"), `[`, 1), sapply(strsplit(zhaoming.etal.vars$KEY.pos,":"), `[`, 2), sapply(strsplit(zhaoming.etal.vars$KEY.pos,":"), `[`, 2))




#############################
#############################
## SNPEFF ##### Annotation ##
#############################
#############################
# https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/")
## read annotated SJLIFE annotated VCF 
# Loop over all chromosomes
chromosomes <- 1:22


## Now extract variants of ClinVar and MetaSVM significance
FINAL.VCF <- {}

capture.output (for( i in 1:length(chromosomes)){
print(paste0("Doing chromosome ", chromosomes[i]))
  
VCF <- fread(paste0("MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr", chromosomes[i], ".PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt"))

#####################################
## Clinvar based P/LP VCF variants ##
#####################################
print(as.data.frame(table(VCF$CLNSIG)))
# Var1     Freq
# 1                                                          . 19388477
# 2                                                    Affects       51
# 3                                                     Benign   117561
# 4                                       Benign/Likely_benign    15102
# 5                                Benign/Likely_benign|_other        7
# 6                                              Benign|_other       14
# 7               Conflicting_interpretations_of_pathogenicity    19244
# 8  Conflicting_interpretations_of_pathogenicity|_association        6
# 9  Conflicting_interpretations_of_pathogenicity|_risk_factor        6
# 10                                             Likely_benign    62840
# 11                                         Likely_pathogenic      593
# 12                                                Pathogenic     2028
# 13                              Pathogenic/Likely_pathogenic      875
# 14                                   Pathogenic|_risk_factor        9
# 15                                    Uncertain_significance    47896
# 16                                               association        7
# 17                                              not_provided      535
# 18                                                protective        3
# 19                                   protective|_risk_factor        5
# 20                                               risk_factor       22

## Wanted Clinvar patterns
WANTED.types.clinvar <- c("^Pathogenic/Likely_pathogenic$|^Likely_pathogenic$|^Likely_pathogenic/Pathogenic$|^Pathogenic$|Pathogenic\\|_risk_factor|^Pathogenic\\|")
print(paste0("Total P or LP vars from clinvar: ", sum(grepl(WANTED.types.clinvar, VCF$CLNSIG, ignore.case = T))))

## Wanted MetaSVM patterns D (available patterns: D= Deleterious; T= Tolerated)
wanted.types.MetaSVM <- c("D")
print(paste0("Total Deleterious vars from MetaSVM: ", sum(grepl(wanted.types.MetaSVM, VCF$dbNSFP_MetaSVM_pred, ignore.case = T))))

VCF.clinvar <- VCF[grepl(WANTED.types.clinvar, VCF$CLNSIG, ignore.case = T),]
VCF.clinvar$PRED_TYPE <- "Clinvar"
VCF.MetaSVM <- VCF[grepl(wanted.types.MetaSVM, VCF$dbNSFP_MetaSVM_pred, ignore.case = T),]
VCF.MetaSVM$PRED_TYPE="MetaSVM"
VCF <- rbind.data.frame(VCF.clinvar, VCF.MetaSVM)
FINAL.VCF <- rbind.data.frame(FINAL.VCF, VCF)
}, file = "SNPEFF_clinvar_metaSVM_from_R_filtering_process.log")

# Qing et al 2020; We considered a missense variant highfunctional impact if
# classified as Deleterious by MetaSVM or listed as Pathogenic/
# Likely-Pathogenic in ClinVar
FINAL.VCF.missense <- FINAL.VCF[grepl("missense", FINAL.VCF$`ANN[*].EFFECT`),]

save.image("SNPEFF_clinvar_metaSVM_from_R_filtering_process.RData")

load("SNPEFF_clinvar_metaSVM_from_R_filtering_process.RData")

#######################################################################
## Check which of the variants from the previous studies are present ##
#######################################################################
# # First remove any leading and trailing spaces
FINAL.VCF$CHROM <- trimws(FINAL.VCF$CHROM, which = "both")
FINAL.VCF$POS <- trimws(FINAL.VCF$POS, which = "both")


## Label indels and SNVs
qin.etal.vars$varTypes <- NULL
qin.etal.vars$varTypes[grepl("^A$|^T$|^G$|^C$", qin.etal.vars$Reference_Allele) & grepl("^A$|^T$|^G$|^C$", qin.etal.vars$Mutant_Allele)] <- "SNV"
qin.etal.vars$CHROM <- paste0("chr",qin.etal.vars$Chr)
qin.etal.vars$START <- as.numeric(qin.etal.vars$Pos_GRCh38)
qin.etal.vars$END <- as.numeric(qin.etal.vars$Pos_GRCh38)
qin.etal.vars$varTypes [!grepl("SNV", qin.etal.vars$varTypes)] <- "INDEL"
qin.etal.vars$END [grepl("INDEL", qin.etal.vars$varTypes)] <- NA
qin.etal.vars$END [grepl("bp", qin.etal.vars$Reference_Allele)] <- as.numeric(gsub("bp| ","", qin.etal.vars$Reference_Allele 
                                                                                   [grepl("bp", qin.etal.vars$Reference_Allele)])) +
  as.numeric(qin.etal.vars$START [grepl("bp", qin.etal.vars$Reference_Allele)]) + 1


qin.etal.vars$END [grepl("bp", qin.etal.vars$Mutant_Allele)] <- as.numeric(gsub("bp| ","", qin.etal.vars$Mutant_Allele 
                                                                                [grepl("bp", qin.etal.vars$Mutant_Allele)])) +
  as.numeric(qin.etal.vars$START [grepl("bp", qin.etal.vars$Mutant_Allele)]) + 1


qin.etal.vars$END[is.na(qin.etal.vars$END)&qin.etal.vars$Reference_Allele !="-"] <- nchar(gsub("-", "", qin.etal.vars$Reference_Allele))[is.na(qin.etal.vars$END)&qin.etal.vars$Reference_Allele !="-"] + 
  + as.numeric(qin.etal.vars$START[is.na(qin.etal.vars$END)&qin.etal.vars$Reference_Allele !="-"]) +1


qin.etal.vars$END[is.na(qin.etal.vars$END)&qin.etal.vars$Mutant_Allele !="-"] <- nchar(gsub("-", "", qin.etal.vars$Mutant_Allele))[is.na(qin.etal.vars$END)&qin.etal.vars$Mutant_Allele !="-"] + 
  + as.numeric(qin.etal.vars$START[is.na(qin.etal.vars$END)&qin.etal.vars$Mutant_Allele !="-"]) +1

qin.etal.vars$START[qin.etal.vars$varTypes == "INDEL"] <- as.numeric(qin.etal.vars$START) [qin.etal.vars$varTypes == "INDEL"]-1

## BED file
write.table(cbind.data.frame(qin.etal.vars$CHROM,qin.etal.vars$START, qin.etal.vars$END), "qin_et_al_variants.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Now check how many of the SNVs are in LOF and CLINVAR annotation in our dataset
qin.SNV <- qin.etal.vars[qin.etal.vars$varTypes == "SNV",]
qin.INDEL <- qin.etal.vars[qin.etal.vars$varTypes == "INDEL",]

length(unique(qin.SNV$KEY.varID))
# 211
sum(unique(qin.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID[FINAL.VCF$PRED_TYPE == "Clinvar"])
# 91
sum(unique(qin.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID[FINAL.VCF$PRED_TYPE == "MetaSVM"])
# 20
## CLINVAR+METASVM
sum(unique(qin.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID)
# 95

# In LOF
LOF.df <- fread("LoF_variants_ID.txt", header = F)
LOF.df$Key.Pos <- sub('^([^:]+:[^:]+).*', '\\1', LOF.df$V1)
LOF <- as.character(LOF.df$V1)
sum(unique(qin.SNV$KEY.varID) %in% LOF)
# 188
# In everything
sum(unique(qin.SNV$KEY.varID) %in% c(LOF, FINAL.VCF$KEY.varID))
# 205




########################################
## Now repeat this for Zhaoming et al ##
########################################
## Label indels and SNVs
zhaoming.etal.vars$varTypes <- NULL
zhaoming.etal.vars$varTypes[grepl("^A$|^T$|^G$|^C$", zhaoming.etal.vars$Reference_Allele) & grepl("^A$|^T$|^G$|^C$", zhaoming.etal.vars$Mutant_Allele)] <- "SNV"
zhaoming.etal.vars$CHROM <- paste0("chr",zhaoming.etal.vars$Chr)
zhaoming.etal.vars$START <- as.numeric(zhaoming.etal.vars$Pos_GRCh38)
zhaoming.etal.vars$END <- as.numeric(zhaoming.etal.vars$Pos_GRCh38)
zhaoming.etal.vars$varTypes [!grepl("SNV", zhaoming.etal.vars$varTypes)] <- "INDEL"
zhaoming.etal.vars$END [grepl("INDEL", zhaoming.etal.vars$varTypes)] <- NA
zhaoming.etal.vars$END [grepl("bp", zhaoming.etal.vars$Reference_Allele)] <- as.numeric(gsub("bp| ","", zhaoming.etal.vars$Reference_Allele 
                                                                                   [grepl("bp", zhaoming.etal.vars$Reference_Allele)])) +
  as.numeric(zhaoming.etal.vars$START [grepl("bp", zhaoming.etal.vars$Reference_Allele)]) + 1


zhaoming.etal.vars$END [grepl("bp", zhaoming.etal.vars$Mutant_Allele)] <- as.numeric(gsub("bp| ","", zhaoming.etal.vars$Mutant_Allele 
                                                                                [grepl("bp", zhaoming.etal.vars$Mutant_Allele)])) +
  as.numeric(zhaoming.etal.vars$START [grepl("bp", zhaoming.etal.vars$Mutant_Allele)]) + 1


zhaoming.etal.vars$END[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Reference_Allele !="-"] <- nchar(gsub("-", "", zhaoming.etal.vars$Reference_Allele))[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Reference_Allele !="-"] + 
  + as.numeric(zhaoming.etal.vars$START[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Reference_Allele !="-"]) +1


zhaoming.etal.vars$END[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Mutant_Allele !="-"] <- nchar(gsub("-", "", zhaoming.etal.vars$Mutant_Allele))[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Mutant_Allele !="-"] + 
  + as.numeric(zhaoming.etal.vars$START[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Mutant_Allele !="-"]) +1

zhaoming.etal.vars$START[zhaoming.etal.vars$varTypes == "INDEL"] <- as.numeric(zhaoming.etal.vars$START) [zhaoming.etal.vars$varTypes == "INDEL"]-1

## BED file
write.table(cbind.data.frame(zhaoming.etal.vars$CHROM,zhaoming.etal.vars$START, zhaoming.etal.vars$END), "zhaoming_et_al_variants.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Now check how many of the SNVs are in LOF and CLINVAR annotation in our dataset
zhaoming.SNV <- zhaoming.etal.vars[zhaoming.etal.vars$varTypes == "SNV",]
zhaoming.INDEL <- zhaoming.etal.vars[zhaoming.etal.vars$varTypes == "INDEL",]

length(unique(zhaoming.SNV$KEY.varID))
# 166
sum(unique(zhaoming.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID[FINAL.VCF$PRED_TYPE == "Clinvar"])
# 114
sum(unique(zhaoming.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID[FINAL.VCF$PRED_TYPE == "MetaSVM"])
# 28
## CLINVAR+METASVM
sum(unique(zhaoming.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID)
# 123

# In LOF
LOF <- fread("LoF_variants_ID.txt")
LOF <- as.character(LOF$V1)
sum(unique(zhaoming.SNV$KEY.varID) %in% LOF)
# 128
# In everything
sum(unique(zhaoming.SNV$KEY.varID) %in% c(LOF, FINAL.VCF$KEY.varID))
# 155

#############################
## NOw checking for INDELS ##
#############################
unique(zhaoming.INDEL$KEY.varID)
zhaoming.INDEL <- zhaoming.INDEL %>% distinct(KEY.varID, .keep_all = TRUE)

zhaoming.INDEL$Pos_GRCh38 <- as.numeric(zhaoming.INDEL$Pos_GRCh38)
zhaoming.INDEL$Pos_GRCh38.minus1 <- zhaoming.INDEL$Pos_GRCh38-1
zhaoming.INDEL$KEY.pos.indels <- paste0("chr", zhaoming.INDEL$Chr, ":", zhaoming.INDEL$Pos_GRCh38.minus1)
write.table(cbind.data.frame(zhaoming.INDEL$CHROM,zhaoming.INDEL$START, zhaoming.INDEL$END), "zhaoming_et_al_variants_INDEL.bed", row.names = F, col.names = F, quote = F, sep = "\t")


unique(qin.INDEL$KEY.varID)
qin.INDEL <- qin.INDEL %>% distinct(KEY.varID, .keep_all = TRUE)

qin.INDEL$Pos_GRCh38 <- as.numeric(qin.INDEL$Pos_GRCh38)
qin.INDEL$Pos_GRCh38.minus1 <- qin.INDEL$Pos_GRCh38-1
qin.INDEL$KEY.pos.indels <- paste0("chr", qin.INDEL$Chr, ":", qin.INDEL$Pos_GRCh38.minus1)
write.table(cbind.data.frame(qin.INDEL$CHROM, qin.INDEL$START, qin.INDEL$END), "qin_et_al_variants_INDEL.bed", row.names = F, col.names = F, quote = F, sep = "\t")










length(unique(zhaoming.INDEL$KEY.pos.indels))
# 130
sum(unique(zhaoming.INDEL$KEY.pos.indels) %in% FINAL.VCF$KEY.pos[FINAL.VCF$PRED_TYPE == "Clinvar"])
# 39
sum(unique(zhaoming.INDEL$KEY.pos.indels) %in% FINAL.VCF$KEY.pos[FINAL.VCF$PRED_TYPE == "MetaSVM"])
# 0
sum(unique(zhaoming.INDEL$KEY.pos.indels) %in% LOF.df$Key.Pos)
# 65
sum(unique(zhaoming.INDEL$KEY.pos.indels) %in% c(FINAL.VCF$KEY.pos, LOF.df$Key.Pos))
# 66



#####################################################################
#####################################################################
#####################################################################
#####################################################################


## Check these in original VCFs from Yadav
zhaoming_original_VCF <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/attr_fraction/zhaoming_vars_vcf.txt")
zhaoming_original_VCF$V2 <- paste0("chr", zhaoming_original_VCF$V1,":", zhaoming_original_VCF$V4, ":", zhaoming_original_VCF$V5, ":", zhaoming_original_VCF$V6)
zhaoming_original_VCF$KEY.Pos <- paste0("chr", zhaoming_original_VCF$V1,":", zhaoming_original_VCF$V4)
write.table(zhaoming_original_VCF, "Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/attr_fraction/zhaoming_vars_vcf_edited.txt", col.names = F, row.names = F, quote = F, sep = "\t")
sum(unique(zhaoming.etal.vars$KEY.pos) %in% zhaoming_original_VCF$KEY.Pos)
sum(unique(zhaoming.etal.vars$KEY.varID) %in% zhaoming_original_VCF$V2)


qin_original_VCF <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/attr_fraction/qin_vars_vcf.txt")
qin_original_VCF$V2 <- paste0("chr", qin_original_VCF$V1,":", qin_original_VCF$V4, ":", qin_original_VCF$V5, ":", qin_original_VCF$V6)
qin_original_VCF$KEY.Pos <- paste0("chr", qin_original_VCF$V1,":", qin_original_VCF$V4)
write.table(qin_original_VCF, "Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/attr_fraction/qin_vars_vcf_edited.txt", col.names = F, row.names = F, quote = F, sep = "\t")
sum(unique(qin.etal.vars$KEY.pos) %in% qin_original_VCF$KEY.Pos)
sum(unique(qin.etal.vars$KEY.varID) %in% qin_original_VCF$V2)
# as.data.frame(table(FINAL.VCF$`ANN[*].EFFECT`))
# 
# 
# CLINVAR <- FINAL.VCF[grepl("Clinvar", FINAL.VCF$PRED_TYPE),]
# MetaSVM <- FINAL.VCF[grepl("MetaSVM", FINAL.VCF$PRED_TYPE),]
# 
# 
# ## Clinvar and MetaSVM
# FINAL.VCF$KEY.pos <- paste(FINAL.VCF$CHROM, FINAL.VCF$POS, sep = ":")
# FINAL.VCF$KEY.varID <- paste(FINAL.VCF$CHROM, FINAL.VCF$POS, FINAL.VCF$REF, FINAL.VCF$ALT, sep = ":")
# 
# 
# ## Zhaoming et al
# length(unique(zhaoming.etal.vars$KEY.pos)) # 295
# sum(unique(zhaoming.etal.vars$KEY.pos) %in% FINAL.VCF$KEY.pos)
# # 126
# 
# length(unique(zhaoming.etal.vars$KEY.varID)) # 299
# sum(unique(zhaoming.etal.vars$KEY.varID) %in% FINAL.VCF$KEY.varID)
# # 123
# 
# ## Qin et al
# length(unique(qin.etal.vars$KEY.pos)) # 389
# sum(unique(qin.etal.vars$KEY.pos) %in% FINAL.VCF$KEY.pos)
# # 94
# length(unique(qin.etal.vars$KEY.varID)) # 392
# sum(unique(qin.etal.vars$KEY.varID) %in% FINAL.VCF$KEY.varID)
# # 95
# 
# ############################
# # check how many samples are from CLINVAR
# ## Zhaoming et al
# length(unique(zhaoming.etal.vars$KEY.pos)) # 295
# sum(unique(zhaoming.etal.vars$KEY.pos) %in% CLINVAR$KEY.pos)
# # 115
# 
# length(unique(zhaoming.etal.vars$KEY.varID)) # 299
# sum(unique(zhaoming.etal.vars$KEY.varID) %in% CLINVAR$KEY.varID)
# # 114
# 
# ## Qin et al
# length(unique(qin.etal.vars$KEY.pos)) # 389
# sum(unique(qin.etal.vars$KEY.pos) %in% CLINVAR$KEY.pos)
# # 91
# length(unique(qin.etal.vars$KEY.varID)) # 392
# sum(unique(qin.etal.vars$KEY.varID) %in% CLINVAR$KEY.varID)
# # 91
# 
# # check how many samples are from MetaSVM
# length(unique(zhaoming.etal.vars$KEY.pos)) # 295
# sum(unique(zhaoming.etal.vars$KEY.pos) %in% MetaSVM$KEY.pos)
# # 30
# 
# length(unique(zhaoming.etal.vars$KEY.varID)) # 299
# sum(unique(zhaoming.etal.vars$KEY.varID) %in% MetaSVM$KEY.varID)
# # 28
# 
# ## Qin et al
# length(unique(qin.etal.vars$KEY.pos)) # 389
# sum(unique(qin.etal.vars$KEY.pos) %in% MetaSVM$KEY.pos)
# # 19
# length(unique(qin.etal.vars$KEY.varID)) # 392
# sum(unique(qin.etal.vars$KEY.varID) %in% MetaSVM$KEY.varID)
# # 20
# 
# write.table(paste0(unique(zhaoming.etal.vars$KEY.varID),":"), "zhaoming_variants.txt", row.names = F, col.names = F, quote = F)
# write.table(paste0(unique(zhaoming.etal.vars$KEY.pos),":"), "zhaoming_variant_sites.txt", row.names = F, col.names = F, quote = F)

