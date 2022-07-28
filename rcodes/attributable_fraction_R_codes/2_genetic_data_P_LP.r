#########################
## Load Phenotype data ##
#########################
# Note: I had to run this on HPC. The datasets are too big to handle

library(data.table)
library(stringr)

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/annotated_variants/annotated_vars_from_PreQC_VCF_Clinvar_MetaSVM_LoF")

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


##############
## ADD P/LP ##
##############
# 1. Variable with carriers status for all Zhaoming and Qin's variants (regardless of P/LP status) have been already added in 1_demographics.r

# 2. Hallmark of Cancer

# 3. All P/LP

#################################
## 2. Hallmark of Cancer genes ##
#################################
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
hallmark.cancer.clinvar.LoF <- rbind.data.frame(hallmark.cancer.clinvar, hallmark.cancer.LoF)
dim(hallmark.cancer.clinvar.LoF)
# 21794
hallmark.cancer.clinvar.LoF <- hallmark.cancer.clinvar.LoF[!duplicated(hallmark.cancer.clinvar.LoF$KEY),]
dim(hallmark.cancer.clinvar.LoF)
# 21367

sum(hallmark.cancer.clinvar.LoF$KEY %in% colnames(P_LP.extracted))
# 21367
hallmark.cancer.clinvar.LoF.PL <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% hallmark.cancer.clinvar.LoF$KEY))]

hallmark.cancer.clinvar.LoF.PL$hallmark.cancer.clinvar.LoF.Non.Ref.Counts <- rowSums(hallmark.cancer.clinvar.LoF.PL[!(colnames(hallmark.cancer.clinvar.LoF.PL) %in% "IID")], na.rm = T) 
hallmark.cancer.clinvar.LoF.PL$hallmark.cancer.clinvar.LoF_carriers <- factor(ifelse(hallmark.cancer.clinvar.LoF.PL$hallmark.cancer.clinvar.LoF.Non.Ref.Counts == 0, "N", "Y"))
hallmark.cancer.clinvar.LoF.PL$hallmark.cancer.clinvar.LoF_carriers  <- factor(hallmark.cancer.clinvar.LoF.PL$hallmark.cancer.clinvar.LoF_carriers, levels = c("N", "Y"))
hallmark.cancer.clinvar.LoF.PL <- cbind.data.frame(IID = hallmark.cancer.clinvar.LoF.PL$IID,
                                                           hallmark.cancer.clinvar.LoF.Non.Ref.Counts = hallmark.cancer.clinvar.LoF.PL$hallmark.cancer.clinvar.LoF.Non.Ref.Counts, 
                                                           hallmark.cancer.clinvar.LoF_carriers = hallmark.cancer.clinvar.LoF.PL$hallmark.cancer.clinvar.LoF_carriers)



# 2. With CliniVar, LoF and MetaSVM
sum(MetaSVM.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity)
# 10311
hallmark.cancer.MetaSVM <- MetaSVM.unique[MetaSVM.unique$`ANN[*].GENE` %in% hallmark.cancer$Adaptive_immunity,]

##
hallmark.cancer.clinvar.LoF.MetaSVM <- rbind.data.frame(hallmark.cancer.clinvar, hallmark.cancer.LoF,hallmark.cancer.MetaSVM)
dim(hallmark.cancer.clinvar.LoF.MetaSVM)
# 32105
hallmark.cancer.clinvar.LoF.MetaSVM <- hallmark.cancer.clinvar.LoF.MetaSVM[!duplicated(hallmark.cancer.clinvar.LoF.MetaSVM$KEY),]
dim(hallmark.cancer.clinvar.LoF.MetaSVM)
# 31117

sum(hallmark.cancer.clinvar.LoF.MetaSVM$KEY %in% colnames(P_LP.extracted))
# 31117
hallmark.cancer.clinvar.LoF.MetaSVM.PL <- P_LP.extracted[c(1, which(colnames(P_LP.extracted) %in% hallmark.cancer.clinvar.LoF.MetaSVM$KEY))]

hallmark.cancer.clinvar.LoF.MetaSVM.PL$hallmark.cancer.clinvar.LoF.MetaSVM.Non.Ref.Counts <- rowSums(hallmark.cancer.clinvar.LoF.MetaSVM.PL[!(colnames(hallmark.cancer.clinvar.LoF.MetaSVM.PL) %in% "IID")], na.rm = T) 
hallmark.cancer.clinvar.LoF.MetaSVM.PL$hallmark.cancer.clinvar.LoF.MetaSVM_carriers <- factor(ifelse(hallmark.cancer.clinvar.LoF.MetaSVM.PL$hallmark.cancer.clinvar.LoF.MetaSVM.Non.Ref.Counts == 0, "N", "Y"))
hallmark.cancer.clinvar.LoF.MetaSVM.PL$hallmark.cancer.clinvar.LoF.MetaSVM_carriers  <- factor(hallmark.cancer.clinvar.LoF.MetaSVM.PL$hallmark.cancer.clinvar.LoF.MetaSVM_carriers, levels = c("N", "Y"))

hallmark.cancer.clinvar.LoF.MetaSVM.PL <- cbind.data.frame(IID = hallmark.cancer.clinvar.LoF.MetaSVM.PL$IID,
                                                           hallmark.cancer.clinvar.LoF.MetaSVM.Non.Ref.Counts = hallmark.cancer.clinvar.LoF.MetaSVM.PL$hallmark.cancer.clinvar.LoF.MetaSVM.Non.Ref.Counts, 
                                                           hallmark.cancer.clinvar.LoF.MetaSVM_carriers = hallmark.cancer.clinvar.LoF.MetaSVM.PL$hallmark.cancer.clinvar.LoF.MetaSVM_carriers)





hallmark.of.cancer.carriers <- cbind.data.frame(hallmark.cancer.clinvar.LoF.PL, hallmark.cancer.clinvar.LoF.MetaSVM.PL[match(hallmark.cancer.clinvar.LoF.PL$IID, hallmark.cancer.clinvar.LoF.MetaSVM.PL$IID),-1])


#####################################
## 3. Pathogenic/Likely Pathogenic ##
#####################################






#########################
## Extract Ethnicities ##
#########################
PHENO.ANY_SN.EUR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]


# lapply(list(wgspop, wgsdiag, subneo, radiation, drug, demog, adultbmi, adolhabits, adlthabits), dim)
# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/demographic.RDATA")
save.image("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/demographic.RDATA")