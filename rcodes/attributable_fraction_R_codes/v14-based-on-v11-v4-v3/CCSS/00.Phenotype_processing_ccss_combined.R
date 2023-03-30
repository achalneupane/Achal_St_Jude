rm()
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_org_Genetic_data_P_LP_v14.Rdata")
ccss_org.PHENO.ANY_SN <- PHENO.ANY_SN
ccss_org.subneo <- subneo

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_exp_Genetic_data_P_LP_v14.Rdata")
ccss_exp.PHENO.ANY_SN <- PHENO.ANY_SN
ccss_exp.subneo <- subneo

sum(colnames(ccss_org.PHENO.ANY_SN) != colnames(ccss_exp.PHENO.ANY_SN))
sum(colnames(ccss_org.subneo) != colnames(ccss_exp.subneo))

PHENO.ANY_SN <- rbind.data.frame(ccss_org.PHENO.ANY_SN, ccss_exp.PHENO.ANY_SN)
PHENO.ANY_SN$ccssid <- sapply(strsplit(PHENO.ANY_SN$ccssid,"_"), `[`, 1)

subneo <- rbind.data.frame(ccss_org.subneo, ccss_exp.subneo)
subneo$ccssid <- sapply(strsplit(subneo$ccssid,"_"), `[`, 1)


## recode categories to make consistently same with SJLIFE categories
# PHENO.ANY_SN$maxsegrtdose.category: None >0-<18 >=18-<30 >=30 Unknown
# PHENO.ANY_SN$maxneckrtdose.category: None >0-<11 >=11-<20 >=20-<30 >=30 Unknown
# PHENO.ANY_SN$maxabdrtdose.category: None >0-<30 >=30 Unknown
# PHENO.ANY_SN$maxchestrtdose.category: None >0-<20 >=20 Unknown
# PHENO.ANY_SN$maxpelvisrtdose.category: None >0-<20 >=20 Unknown


PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
# None 0-18 18-30 >=30 Unknown -->>> None >0-<18 >=18-<30 >=30 Unknown
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "0-18"] <- ">0-<18"
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "18-30"] <- ">=18-<30"
PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("None", ">0-<18", ">=18-<30", ">=30", "Unknown"))
table(PHENO.ANY_SN$maxsegrtdose.category)

PHENO.ANY_SN$maxneckrtdose.category <- as.character(PHENO.ANY_SN$maxneckrtdose.category)
# None 0-11 11-20 20-30 >=30 Unknown -->>> None >0-<11 >=11-<20 >=20-<30 >=30 Unknown
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "0-11"] <- ">0-<11"
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "11-20"] <- ">=11-<20"
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "20-30"] <- ">=20-<30"
PHENO.ANY_SN$maxneckrtdose.category <- factor(PHENO.ANY_SN$maxneckrtdose.category, levels = c("None", ">0-<11", ">=11-<20", ">=20-<30", ">=30", "Unknown"))
table(PHENO.ANY_SN$maxneckrtdose.category)

PHENO.ANY_SN$maxabdrtdose.category <- as.character(PHENO.ANY_SN$maxabdrtdose.category)
#  None 0-30 >=30 Unknown -->>> None >0-<30 >=30 Unknown
PHENO.ANY_SN$maxabdrtdose.category[PHENO.ANY_SN$maxabdrtdose.category == "0-30"] <- ">0-<30"
PHENO.ANY_SN$maxabdrtdose.category <- factor(PHENO.ANY_SN$maxabdrtdose.category, levels = c("None", ">0-<30", ">=30", "Unknown"))
table(PHENO.ANY_SN$maxabdrtdose.category)

PHENO.ANY_SN$maxchestrtdose.category <- as.character(PHENO.ANY_SN$maxchestrtdose.category)
#  None 0-20 >=20 Unknown -->>> None >0-<20 >=20 Unknown
PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == "0-20"] <- ">0-<20"
PHENO.ANY_SN$maxchestrtdose.category <- factor(PHENO.ANY_SN$maxchestrtdose.category, levels = c("None", ">0-<20", ">=20", "Unknown"))
table(PHENO.ANY_SN$maxchestrtdose.category)

PHENO.ANY_SN$maxpelvisrtdose.category <- as.character(PHENO.ANY_SN$maxpelvisrtdose.category)
#  None 0-20 >=20 Unknown -->>> None >0-<20 >=20 Unknown
PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == "0-20"] <- ">0-<20"
PHENO.ANY_SN$maxpelvisrtdose.category <- factor(PHENO.ANY_SN$maxpelvisrtdose.category, levels = c("None", ">0-<20", ">=20", "Unknown"))
table(PHENO.ANY_SN$maxpelvisrtdose.category)

rm(list=ls()[!grepl(c("PHENO.ANY_SN|subneo"), ls())])

save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_combined_Genetic_data_P_LP_v14.Rdata")
