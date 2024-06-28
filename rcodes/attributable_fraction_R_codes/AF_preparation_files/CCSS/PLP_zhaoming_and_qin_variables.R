rm(list=ls())

#########
## PLP ##
#########
## Read sample map file
mapfile <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt", header = F, sep = " ", stringsAsFactors = F)

## Zhaoming 
Sys.setlocale("LC_ALL", "C")
zhaoming_tables <- read.delim("/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Zhaoming_Wang_Genetic_risk_for_subsequent_neoplasms_JCO.2018.77.8589/P-PL-sjlife-genetics-sn_SNV_INDELS_07_13_2022.txt", sep = "\t")
dim(zhaoming_tables)
sum(!duplicated(zhaoming_tables$FULLKEY))

## Read file with corresponding CCSS variants
zhaoming.ccss.var.id <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/Zhaoming_Wang_et_al_variants.txt", header = T, sep = "\t")
table(zhaoming_tables$FULLKEY %in% zhaoming.ccss.var.id$KEY.varID)
zhaoming_tables$varID.our.data <- zhaoming.ccss.var.id$ID[match(zhaoming_tables$FULLKEY, zhaoming.ccss.var.id$KEY.varID)]

# Keeping variants in Table A4 only
tableA4 <- zhaoming_tables[grepl ("A4", zhaoming_tables$Table),]
dim(tableA4)
# 166

## Read raw file from ccss Zhaoming
zhaoming.ccss.raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/zhaoming_ccss/merged_data_plink_recodeA.raw", header = T)
zhaoming.ccss.raw <- zhaoming.ccss.raw[-grep("FID|PAT|MAT|SEX|PHENOTYPE", colnames(zhaoming.ccss.raw))]
# edit column names
colnames(zhaoming.ccss.raw) <- str_split(gsub("\\.", ":", colnames(zhaoming.ccss.raw)), "_", simplify=T)[,1]
table(colnames(zhaoming.ccss.raw) %in% tableA4$varID.our.data)

## Keep CCSS_exp samples only
zhaoming.ccss.raw$IID <- mapfile$V4[match(zhaoming.ccss.raw$IID, mapfile$V3)]
zhaoming.ccss.raw <- zhaoming.ccss.raw[!(grepl("SJL", zhaoming.ccss.raw$IID) | is.na(zhaoming.ccss.raw$IID)),]

# Keeping variants in Table A4 only
zhaoming.ccss.raw <- zhaoming.ccss.raw[, c(1, which(colnames(zhaoming.ccss.raw) %in% tableA4$varID.our.data))]
zhaoming.ccss.raw$Zhaoming_Non.Ref.Counts <- rowSums(zhaoming.ccss.raw[-1]) 
zhaoming.ccss.raw$Zhaoming_carriers <- factor(ifelse(zhaoming.ccss.raw$Zhaoming_Non.Ref.Counts == 0, "N", "Y"))

# pp <- search_list

## Qin
qin.var.id <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/Na_Qin_et_al_variants.txt", header = T, sep = "\t")

qin.ccss.raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/genetic_data/qin_ccss/merged_data_plink_recodeA.raw", header = T)
qin.ccss.raw <- qin.ccss.raw[-grep("FID|PAT|MAT|SEX|PHENOTYPE", colnames(qin.ccss.raw))]
colnames(qin.ccss.raw) <- str_split(gsub("\\.", ":", colnames(qin.ccss.raw)), "_", simplify=T)[,1]
dim(qin.ccss.raw)
# 7979  416

## Keep CCSS_exp samples only
qin.ccss.raw$IID <- mapfile$V4[match(qin.ccss.raw$IID, mapfile$V3)]
qin.ccss.raw <- qin.ccss.raw[!(grepl("SJL", qin.ccss.raw$IID) | is.na(qin.ccss.raw$IID)),]



sum(colnames(qin.ccss.raw) %in% qin.var.id$ID )
# 386
qin.ccss.raw <- qin.ccss.raw[, c(1, which(colnames(qin.ccss.raw) %in% qin.var.id$ID))]

qin.ccss.raw <- qin.ccss.raw[,c(1,which(!colnames(qin.ccss.raw) %in% colnames(zhaoming.ccss.raw)))]
dim(qin.ccss.raw)
# 3007  349

qin.ccss.raw$Qin_without_Zhaoming_vars.Non.Ref.Counts <- rowSums(qin.ccss.raw[-1]) 
qin.ccss.raw$Qin_without_Zhaoming_vars_carriers <- factor(ifelse(qin.ccss.raw$Qin_without_Zhaoming_vars.Non.Ref.Counts == 0, "N", "Y"))


# Keep CCSS expansion only
ccss_exp.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/merged.dat.fam")

ccss_exp.samples$Zhaoming_carriers <-  zhaoming.ccss.raw$Zhaoming_carriers[match(ccss_exp.samples$V2, zhaoming.ccss.raw$IID)]
ccss_exp.samples$Qin_without_Zhaoming_vars_carriers <-  qin.ccss.raw$Qin_without_Zhaoming_vars_carriers[match(ccss_exp.samples$V2, qin.ccss.raw$IID)]

table(ccss_exp.samples$Zhaoming_carriers)
# 16
table(ccss_exp.samples$Qin_without_Zhaoming_vars_carriers)
# 94

save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss_exp_P_LP_zhaoming_qin_v17.Rdata")

