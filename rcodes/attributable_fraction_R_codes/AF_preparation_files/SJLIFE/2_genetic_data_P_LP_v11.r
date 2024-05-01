#########################
## Load Phenotype data ##
#########################
# Note: I had to run this on HPC. The datasets are too big to handle; Some of the variables loaded here are derived from this code: variants_counts_preQC_VCF_06_13_2022.R

library(data.table)
library(stringr)

# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/annotated_variants/annotated_vars_from_PreQC_VCF_Clinvar_MetaSVM_LoF")

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/1_demographics.RDATA")
load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/1_demographics.RDATA")

# Also load P/LP data from pre-QC VCF
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")
load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")

# read all P/LP extracted from VCF. These need to be filtered by checking the alleles
# P_LP.extracted <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/annotated_variants/annotated_vars_from_PreQC_VCF_Clinvar_MetaSVM_LoF/sjlife_all_pathogenic_variants_PreQC_recodeA.raw", header = T)
P_LP.extracted <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/annotated_variants/annotated_vars_from_PreQC_VCF_Clinvar_MetaSVM_LoF/sjlife_all_pathogenic_variants_PreQC_recodeA.raw", header = T)
P_LP.extracted <- as.data.frame(P_LP.extracted)
P_LP.extracted <- P_LP.extracted[-grep("FID|PAT|MAT|SEX|PHENOTYPE", colnames(P_LP.extracted))]
# edit column names
colnames(P_LP.extracted) <- str_split(gsub("\\.", ":", colnames(P_LP.extracted)), "_", simplify=T)[,1]

################################
## Zhaoming and Qins variants ##
################################
zhaoming.qin.variants <- unique(c(colnames(Zhaoming_vars), colnames(QIN_vars)))
zhaoming.qin.variants <- zhaoming.qin.variants[grepl("^chr", zhaoming.qin.variants)]

##############
## ADD P/LP ##
##############
# 1. Variable with carriers status for all Zhaoming and Qin's variants (regardless of P/LP status) have been already added in 1_demographics.r

# 2. Hallmark of Cancer

# 3. All P/LP

# Abbreviations: 
# hallmark.cancer.clinvar.LoF=H.C.Clin.LoF
# without.Zhaoming.Qin = WO.Zhao.Qin
##################################
## 2.a Hallmark of Cancer genes ##
##################################
## 1. With CliniVar and LOF
# hallmark.cancer <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/41467_2020_16293_MOESM4_ESM.txt", sep = "\t")
hallmark.cancer <- read.delim("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/41467_2020_16293_MOESM4_ESM.txt", sep = "\t")
dim(hallmark.cancer)
# 1558
sum(CLINVAR.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity)
# 622
hallmark.cancer.clinvar <- CLINVAR.unique[CLINVAR.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity,]

sum(LoF.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity)
# 21172
hallmark.cancer.LoF <- LoF.unique[LoF.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity,]
H.C.Clin.LoF <- rbind.data.frame(hallmark.cancer.clinvar, hallmark.cancer.LoF)
dim(H.C.Clin.LoF)
# 21794
H.C.Clin.LoF <- H.C.Clin.LoF[!duplicated(H.C.Clin.LoF$KEY),]
dim(H.C.Clin.LoF)
# 21367

sum(H.C.Clin.LoF$KEY %in% colnames(P_LP.extracted))
# 21367
H.C.Clin.LoF.PL <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% H.C.Clin.LoF$KEY))]

H.C.Clin.LoF.PL$H.C.Clin.LoF.Non.Ref.Counts <- rowSums(H.C.Clin.LoF.PL[!(colnames(H.C.Clin.LoF.PL) %in% "IID")], na.rm = T) 
H.C.Clin.LoF.PL$H.C.Clin.LoF_carriers <- factor(ifelse(H.C.Clin.LoF.PL$H.C.Clin.LoF.Non.Ref.Counts == 0, "N", "Y"))
H.C.Clin.LoF.PL$H.C.Clin.LoF_carriers  <- factor(H.C.Clin.LoF.PL$H.C.Clin.LoF_carriers, levels = c("N", "Y"))
H.C.Clin.LoF.PL <- cbind.data.frame(IID = H.C.Clin.LoF.PL$IID,
                        H.C.Clin.LoF.Non.Ref.Counts = H.C.Clin.LoF.PL$H.C.Clin.LoF.Non.Ref.Counts, 
                        H.C.Clin.LoF_carriers = H.C.Clin.LoF.PL$H.C.Clin.LoF_carriers)



# 2. With CliniVar, LoF and MetaSVM
sum(MetaSVM.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity)
# 10311
hallmark.cancer.MetaSVM <- MetaSVM.unique[MetaSVM.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity,]

##
H.C.Clin.LoF.MetaSVM <- rbind.data.frame(hallmark.cancer.clinvar, hallmark.cancer.LoF,hallmark.cancer.MetaSVM)
dim(H.C.Clin.LoF.MetaSVM)
# 32105
H.C.Clin.LoF.MetaSVM <- H.C.Clin.LoF.MetaSVM[!duplicated(H.C.Clin.LoF.MetaSVM$KEY),]
dim(H.C.Clin.LoF.MetaSVM)
# 31117

sum(H.C.Clin.LoF.MetaSVM$KEY %in% colnames(P_LP.extracted))
# 31117
H.C.Clin.LoF.MetaSVM.PL <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% H.C.Clin.LoF.MetaSVM$KEY))]

H.C.Clin.LoF.MetaSVM.PL$H.C.Clin.LoF.MetaSVM.Non.Ref.Counts <- rowSums(H.C.Clin.LoF.MetaSVM.PL[!(colnames(H.C.Clin.LoF.MetaSVM.PL) %in% "IID")], na.rm = T) 
H.C.Clin.LoF.MetaSVM.PL$H.C.Clin.LoF.MetaSVM_carriers <- factor(ifelse(H.C.Clin.LoF.MetaSVM.PL$H.C.Clin.LoF.MetaSVM.Non.Ref.Counts == 0, "N", "Y"))
H.C.Clin.LoF.MetaSVM.PL$H.C.Clin.LoF.MetaSVM_carriers  <- factor(H.C.Clin.LoF.MetaSVM.PL$H.C.Clin.LoF.MetaSVM_carriers, levels = c("N", "Y"))

H.C.Clin.LoF.MetaSVM.PL <- cbind.data.frame(IID = H.C.Clin.LoF.MetaSVM.PL$IID,
                                H.C.Clin.LoF.MetaSVM.Non.Ref.Counts = H.C.Clin.LoF.MetaSVM.PL$H.C.Clin.LoF.MetaSVM.Non.Ref.Counts, 
                                H.C.Clin.LoF.MetaSVM_carriers = H.C.Clin.LoF.MetaSVM.PL$H.C.Clin.LoF.MetaSVM_carriers)


hallmark.of.cancer.carriers <- cbind.data.frame(H.C.Clin.LoF.PL, H.C.Clin.LoF.MetaSVM.PL[match(H.C.Clin.LoF.PL$IID, H.C.Clin.LoF.MetaSVM.PL$IID),-1])


###########################################
## 2.b Without zhaoming and Qin variants ##
###########################################
## 1. With CliniVar and LOF
# hallmark.cancer <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/41467_2020_16293_MOESM4_ESM.txt", sep = "\t")
hallmark.cancer <- read.delim("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/41467_2020_16293_MOESM4_ESM.txt", sep = "\t")
dim(hallmark.cancer)
# 1558
sum(CLINVAR.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity)
# 622
hallmark.cancer.clinvar <- CLINVAR.unique[CLINVAR.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity,]

sum(LoF.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity)
# 21172
hallmark.cancer.LoF <- LoF.unique[LoF.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity,]
H.C.Clin.LoF <- rbind.data.frame(hallmark.cancer.clinvar, hallmark.cancer.LoF)
dim(H.C.Clin.LoF)
# 21794
H.C.Clin.LoF <- H.C.Clin.LoF[!duplicated(H.C.Clin.LoF$KEY),]
dim(H.C.Clin.LoF)
# 21367

sum(H.C.Clin.LoF$KEY %in% colnames(P_LP.extracted))
# 21367


H.C.Clin.LoF.PL <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% H.C.Clin.LoF$KEY))]

## Exclude Zhaoming's and Qin's variants
H.C.Clin.LoF.PL.WO.Zhao.Qin.variants <- H.C.Clin.LoF.PL[!colnames(H.C.Clin.LoF.PL) %in% zhaoming.qin.variants]

H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts <- rowSums(H.C.Clin.LoF.PL.WO.Zhao.Qin.variants[!(colnames(H.C.Clin.LoF.PL.WO.Zhao.Qin.variants) %in% "IID")], na.rm = T) 
H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.WO.Zhao.Qin.variants_carriers <- factor(ifelse(H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts == 0, "N", "Y"))
H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.WO.Zhao.Qin.variants_carriers  <- factor(H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.WO.Zhao.Qin.variants_carriers, levels = c("N", "Y"))
H.C.Clin.LoF.PL.WO.Zhao.Qin.variants <- cbind.data.frame(IID = H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$IID,
                                             H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts = H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.WO.Zhao.Qin_variants.Non.Ref.Counts, 
                                             H.C.Clin.LoF.WO.Zhao.Qin.variants_carriers = H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.WO.Zhao.Qin.variants_carriers)



# 2. With CliniVar, LoF and MetaSVM
sum(MetaSVM.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity)
# 10311
hallmark.cancer.MetaSVM <- MetaSVM.unique[MetaSVM.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity,]

##
H.C.Clin.LoF.MetaSVM <- rbind.data.frame(hallmark.cancer.clinvar, hallmark.cancer.LoF,hallmark.cancer.MetaSVM)
dim(H.C.Clin.LoF.MetaSVM)
# 32105
H.C.Clin.LoF.MetaSVM <- H.C.Clin.LoF.MetaSVM[!duplicated(H.C.Clin.LoF.MetaSVM$KEY),]
dim(H.C.Clin.LoF.MetaSVM)
# 31117

sum(H.C.Clin.LoF.MetaSVM$KEY %in% colnames(P_LP.extracted))
# 31117
H.C.Clin.LoF.MetaSVM.PL <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% H.C.Clin.LoF.MetaSVM$KEY))]

## Exclude Zhaoming's and Qin's variants
H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants <- H.C.Clin.LoF.MetaSVM.PL[!colnames(H.C.Clin.LoF.MetaSVM.PL) %in% zhaoming.qin.variants]

H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts <- rowSums(H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants[!(colnames(H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants) %in% "IID")], na.rm = T) 
H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin.variants_carriers <- factor(ifelse(H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts == 0, "N", "Y"))
H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin.variants_carriers  <- factor(H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin.variants_carriers, levels = c("N", "Y"))
H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants <- cbind.data.frame(IID = H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$IID,
                                                     H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts = H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin_variants.Non.Ref.Counts, 
                                                     H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin.variants_carriers = H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$H.C.Clin.LoF.MetaSVM.WO.Zhao.Qin.variants_carriers)





hallmark.of.cancer_WO.Zhao.Qin.variants.carriers <- cbind.data.frame(H.C.Clin.LoF.PL.WO.Zhao.Qin.variants, H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants[match(H.C.Clin.LoF.PL.WO.Zhao.Qin.variants$IID, H.C.Clin.LoF.MetaSVM.PL.WO.Zhao.Qin.variants$IID),-1])

hallmark.of.cancer.carriers <- cbind.data.frame(hallmark.of.cancer.carriers, hallmark.of.cancer_WO.Zhao.Qin.variants.carriers[match(hallmark.of.cancer.carriers$IID, hallmark.of.cancer_WO.Zhao.Qin.variants.carriers$IID),-1])



######################################
## 3.a Pathogenic/Likely Pathogenic ##
######################################
# 1. Clinvar and LoF
All.P.LP.clinvars.LoF <- unique(c(CLINVAR.unique$KEY, LoF.unique$KEY))
sum(All.P.LP.clinvars.LoF %in% colnames(P_LP.extracted))
# 318253
All.P.LP.clinvars.LoF <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% All.P.LP.clinvars.LoF))]
All.P.LP.clinvars.LoF$All.P.LP.clinvars.LoF.Non.Ref.Counts <- rowSums(All.P.LP.clinvars.LoF[!(colnames(All.P.LP.clinvars.LoF) %in% "IID")], na.rm = T) 
All.P.LP.clinvars.LoF$All.P.LP.clinvars.LoF_carriers <- factor(ifelse(All.P.LP.clinvars.LoF$All.P.LP.clinvars.LoF.Non.Ref.Counts == 0, "N", "Y"))
All.P.LP.clinvars.LoF$All.P.LP.clinvars.LoF_carriers  <- factor(All.P.LP.clinvars.LoF$All.P.LP.clinvars.LoF_carriers, levels = c("N", "Y"))
All.P.LP.clinvars.LoF <- cbind.data.frame(IID = All.P.LP.clinvars.LoF$IID, All.P.LP.clinvars.LoF.Non.Ref.Counts = All.P.LP.clinvars.LoF$All.P.LP.clinvars.LoF.Non.Ref.Counts, All.P.LP.clin.LoF_carriers = All.P.LP.clinvars.LoF$All.P.LP.clinvars.LoF_carriers)


# 2. Clinvar, LoF and MetaSVM                          
All.P.LP.clinvars.LoF.MetaSVM <- unique(c(CLINVAR.unique$KEY, LoF.unique$KEY, MetaSVM.unique$KEY))
sum(All.P.LP.clinvars.LoF.MetaSVM %in% colnames(P_LP.extracted))
# 397036
All.P.LP.clinvars.LoF.MetaSVM <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% All.P.LP.clinvars.LoF.MetaSVM))]
All.P.LP.clinvars.LoF.MetaSVM$All.P.LP.clinvars.LoF.MetaSVM.Non.Ref.Counts <- rowSums(All.P.LP.clinvars.LoF.MetaSVM[!(colnames(All.P.LP.clinvars.LoF.MetaSVM) %in% "IID")], na.rm = T) 
All.P.LP.clinvars.LoF.MetaSVM$All.P.LP.clinvars.LoF.MetaSVM_carriers <- factor(ifelse(All.P.LP.clinvars.LoF.MetaSVM$All.P.LP.clinvars.LoF.MetaSVM.Non.Ref.Counts == 0, "N", "Y"))
All.P.LP.clinvars.LoF.MetaSVM$All.P.LP.clinvars.LoF.MetaSVM_carriers  <- factor(All.P.LP.clinvars.LoF.MetaSVM$All.P.LP.clinvars.LoF.MetaSVM_carriers, levels = c("N", "Y"))
All.P.LP.clinvars.LoF.MetaSVM  <- cbind.data.frame(IID = All.P.LP.clinvars.LoF.MetaSVM$IID, All.P.LP.clinvars.LoF.MetaSVM.Non.Ref.Counts = All.P.LP.clinvars.LoF.MetaSVM$All.P.LP.clinvars.LoF.MetaSVM.Non.Ref.Counts, All.P.LP.clin.LoF.MetaSVM_carriers = All.P.LP.clinvars.LoF.MetaSVM$All.P.LP.clinvars.LoF.MetaSVM_carriers)

All.P.LP <- cbind.data.frame(All.P.LP.clinvars.LoF, All.P.LP.clinvars.LoF.MetaSVM[match(All.P.LP.clinvars.LoF$IID, All.P.LP.clinvars.LoF.MetaSVM$IID),-1])

####################################################################
## 3.b Pathogenic/Likely Pathogenic without any previous variants ##
####################################################################
# 1. Clinvar and LoF
All.P.LP.clinvars.LoF.WO.Prior.vars <- unique(c(CLINVAR.unique$KEY, LoF.unique$KEY))
All.P.LP.clinvars.LoF.WO.Prior.vars <- All.P.LP.clinvars.LoF.WO.Prior.vars[!All.P.LP.clinvars.LoF.WO.Prior.vars %in% unique(c(colnames(H.C.Clin.LoF.PL), zhaoming.qin.variants))]
sum(All.P.LP.clinvars.LoF.WO.Prior.vars %in% colnames(P_LP.extracted))
# 296620

All.P.LP.clinvars.LoF.WO.Prior.vars <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% All.P.LP.clinvars.LoF.WO.Prior.vars))]
All.P.LP.clinvars.LoF.WO.Prior.vars$All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts <- rowSums(All.P.LP.clinvars.LoF.WO.Prior.vars[!(colnames(All.P.LP.clinvars.LoF.WO.Prior.vars) %in% "IID")], na.rm = T) 
All.P.LP.clinvars.LoF.WO.Prior.vars$All.P.LP.clinvars.LoF.WO.Prior.vars_carriers <- factor(ifelse(All.P.LP.clinvars.LoF.WO.Prior.vars$All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts == 0, "N", "Y"))
All.P.LP.clinvars.LoF.WO.Prior.vars$All.P.LP.clinvars.LoF.WO.Prior.vars_carriers  <- factor(All.P.LP.clinvars.LoF.WO.Prior.vars$All.P.LP.clinvars.LoF.WO.Prior.vars_carriers, levels = c("N", "Y"))
All.P.LP.clinvars.LoF.WO.Prior.vars <- cbind.data.frame(IID = All.P.LP.clinvars.LoF.WO.Prior.vars$IID, All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts = All.P.LP.clinvars.LoF.WO.Prior.vars$All.P.LP.clinvars.LoF.WO.Prior.vars.Non.Ref.Counts, All.P.LP.clinvars.LoF.WO.Prior.vars_carriers = All.P.LP.clinvars.LoF.WO.Prior.vars$All.P.LP.clinvars.LoF.WO.Prior.vars_carriers)





# 2. Clinvar, LoF and MetaSVM        
All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars <- unique(c(CLINVAR.unique$KEY, LoF.unique$KEY, MetaSVM.unique$KEY))
All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars <- All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars[!All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars %in% unique(c(colnames(H.C.Clin.LoF.MetaSVM.PL), zhaoming.qin.variants))]
sum(All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars %in% colnames(P_LP.extracted))
# 365650

All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars))]
All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts <- rowSums(All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars[!(colnames(All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars) %in% "IID")], na.rm = T) 
All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars_carriers <- factor(ifelse(All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts == 0, "N", "Y"))
All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars_carriers  <- factor(All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars_carriers, levels = c("N", "Y"))
All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars <- cbind.data.frame(IID = All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$IID, All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts = All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars.Non.Ref.Counts, All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars_carriers = All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars_carriers)


All.P.LP.WO.Prior.vars <- cbind.data.frame(All.P.LP.clinvars.LoF.WO.Prior.vars, All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars[match(All.P.LP.clinvars.LoF.WO.Prior.vars$IID, All.P.LP.clinvars.LoF.MetaSVM.WO.Prior.vars$IID),-1])


All.P.LP.carriers <- cbind.data.frame(All.P.LP, All.P.LP.WO.Prior.vars[match(All.P.LP$IID, All.P.LP.WO.Prior.vars$IID), -1])

final.genetic.data <- cbind.data.frame(hallmark.of.cancer.carriers, All.P.LP.carriers[match(hallmark.of.cancer.carriers$IID, All.P.LP.carriers$IID), -1])


# Merge this to Phenotype
PHENO.ANY_SN.with.genetic.data <- cbind.data.frame(PHENO.ANY_SN, final.genetic.data[match(PHENO.ANY_SN$sjlid, final.genetic.data$IID), -1])

## Save RDS
# saveRDS(PHENO.ANY_SN.with.genetic.data, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/2_genetic_data_P_LP.rds")
