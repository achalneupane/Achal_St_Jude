# rm(list=ls())
# 
# ## CCSS samples with WGS; overlapping samples in SJLIFE removed
# length(ccss_samples)
# # 7943
# 
# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2")
# 
# ## Read in CCSS Phenotype/diagnosis from Huiqi
# CCSS_data.ALL <- read.delim("ExportedCCSS_data.txt", header = T, sep = "\t", stringsAsFactors = F)
# ## Keep CCSS with WGS only
# CCSS_data <- CCSS_data.ALL[CCSS_data.ALL$ccssid %in% ccss_samples,]
# nrow(CCSS_data)
# # 9747
# 
# ## Read in Lifestyle from Huiqi
# BMI.PA.SMK.DRK.ALL <- read.delim("ExportedCCSS_BMI_PA_Smk_drink.txt", header = T, sep = "\t", stringsAsFactors = F)
# 
# BMI.PA.SMK.DRK <- BMI.PA.SMK.DRK.ALL[c("ccssid", "a_base", "a_fu1", "a_fu2", "a_fu3", "a_fu2007", "a_fu5", "a_fu6",
#                           "cbmi_0",  "cbmi_2", "cbmi_2007", "cbmi_5",
#                           "t_eqmodfu2", "cdc_fu2", "t_eqmodfu5", "cdc_fu5", "t_eqmodfu6", "cdc_fu6",
#                           "smkcatb", "smkcatf2", "smkcatf07", "smkcatf5",
#                           "riskyb", "riskyf7", "riskyfu5")]
# 
# colnames(BMI.PA.SMK.DRK) <- c("ccssid", "base.age", "fu1.age", "fu2.age", "fu3.age", "fu7.age", "fu5.age", "fu6.age",
#                           "base.bmi",  "fu2.bmi", "fu7.bmi", "fu5.bmi",
#                           "fu2.MET", "fu2.CDC", "fu5.MET", "fu5.CDC", "fu6.MET", "fu6.CDC",
#                           "base.smk", "fu2.smk", "fu7.smk", "fu5.smk",
#                           "base.riskydrk", "fu7.riskydrk", "fu5.riskydrk")



save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2/CCSS_data.Rdata")


# Ensure consistency across follow-ups by adding missing columns
BMI.PA.SMK.DRK$base.MET <- NA
BMI.PA.SMK.DRK$base.CDC <- NA

BMI.PA.SMK.DRK$fu1.bmi <- NA
BMI.PA.SMK.DRK$fu1.MET <- NA
BMI.PA.SMK.DRK$fu1.CDC <- NA
BMI.PA.SMK.DRK$fu1.smk <- NA
BMI.PA.SMK.DRK$fu1.riskydrk <- NA

BMI.PA.SMK.DRK$fu2.riskydrk <- NA

BMI.PA.SMK.DRK$fu3.bmi <- NA
BMI.PA.SMK.DRK$fu3.bmi <- NA
BMI.PA.SMK.DRK$fu3.MET <- NA
BMI.PA.SMK.DRK$fu3.CDC <- NA
BMI.PA.SMK.DRK$fu3.smk <- NA
BMI.PA.SMK.DRK$fu3.riskydrk <- NA

BMI.PA.SMK.DRK$fu7.MET <- NA
BMI.PA.SMK.DRK$fu7.CDC <- NA

BMI.PA.SMK.DRK$fu6.bmi <- NA
BMI.PA.SMK.DRK$fu6.smk <- NA
BMI.PA.SMK.DRK$fu6.riskydrk <- NA

## Replace "." with NA
BMI.PA.SMK.DRK[BMI.PA.SMK.DRK == "."] <- NA


## make prefixes to suffixes
names(BMI.PA.SMK.DRK) <- strsplit(names(BMI.PA.SMK.DRK), '\\.') |> lapply(rev) |> sapply(paste, collapse='.')

# Reshape to long format
ALL.LIFESTYLE <- reshape(BMI.PA.SMK.DRK, direction='l', idvar='ccssid', varying=sort(names(BMI.PA.SMK.DRK)[-1]))

ALL.LIFESTYLE$age <- as.numeric(ALL.LIFESTYLE$age)

## remove rows that have all 4 lifestyle missing
ALL.LIFESTYLE <- ALL.LIFESTYLE %>%
  filter(!(is.na(bmi) & is.na(CDC) & is.na(riskydrk) & is.na(smk)))

## Keep adult only
ALL.LIFESTYLE <- subset(ALL.LIFESTYLE, age >= 18)


# Keep only those that have CCSS WGS
ALL.LIFESTYLE <- ALL.LIFESTYLE[ALL.LIFESTYLE$ccssid %in% ccss_samples,]

######################
## recode lifestyle ##
######################
# Obesity
ALL.LIFESTYLE$Obese_yn <- ifelse(as.numeric(ALL.LIFESTYLE$bmi) < 30, "No", "Yes")
ALL.LIFESTYLE$Obese_yn[is.na(ALL.LIFESTYLE$Obese_yn)] <- "Unknown"
ALL.LIFESTYLE$Obese_yn <- factor(ALL.LIFESTYLE$Obese_yn, level = c("No", "Yes", "Unknown")) 


# Physical activity
ALL.LIFESTYLE$PhysicalActivity_yn <- ifelse(ALL.LIFESTYLE$CDC == "Yes", "Yes", "No")
ALL.LIFESTYLE$PhysicalActivity_yn[is.na(ALL.LIFESTYLE$PhysicalActivity_yn)] <- "Unknown"
ALL.LIFESTYLE$PhysicalActivity_yn <- factor(ALL.LIFESTYLE$PhysicalActivity_yn, level = c("Yes", "No", "Unknown")) 


# Smoker
ALL.LIFESTYLE$Smoker_ever_yn <- ifelse(ALL.LIFESTYLE$smk != 3, "No", "Yes")
ALL.LIFESTYLE$Smoker_ever_yn[is.na(ALL.LIFESTYLE$Smoker_ever_yn)] <- "Unknown"
ALL.LIFESTYLE$Smoker_ever_yn <- factor(ALL.LIFESTYLE$Smoker_ever_yn, level = c("No", "Yes", "Unknown")) 

# risky or heavy drinker
ALL.LIFESTYLE$RiskyHeavyDrink_yn <- ALL.LIFESTYLE$riskydrk 
ALL.LIFESTYLE$RiskyHeavyDrink_yn[is.na(ALL.LIFESTYLE$RiskyHeavyDrink_yn)] <- "Unknown"
ALL.LIFESTYLE$RiskyHeavyDrink_yn <- factor(ALL.LIFESTYLE$RiskyHeavyDrink_yn, level = c("No", "Yes", "Unknown")) 

## Lifestyle variables
nrow(ALL.LIFESTYLE)
# 21115
###########################
## relabel age variables ##
###########################
## Use this variable CCSS_data_ANY_SN for ANY_SN
CCSS_data$AGE.ANY_SN <- as.numeric(CCSS_data$a_candx)
CCSS_data$agedx <- as.numeric(CCSS_data$a_dx)
CCSS_data$agelstcontact <- as.numeric(CCSS_data$a_end)


## Age at diagnosis
# PHENO.ANY_SN$agedx <- floor(PHENO.ANY_SN$agedx)
CCSS_data$AGE_AT_DIAGNOSIS[CCSS_data$agedx >= 0 & CCSS_data$agedx < 5 ] <- "0-4"
CCSS_data$AGE_AT_DIAGNOSIS[CCSS_data$agedx >= 5 & CCSS_data$agedx < 10 ] <- "5-9"
CCSS_data$AGE_AT_DIAGNOSIS[CCSS_data$agedx >= 10 & CCSS_data$agedx < 15 ] <- "10-14"
CCSS_data$AGE_AT_DIAGNOSIS[CCSS_data$agedx >= 15 ] <- ">=15"
CCSS_data$AGE_AT_DIAGNOSIS <- factor(CCSS_data$AGE_AT_DIAGNOSIS, levels = c("0-4", "5-9", "10-14", ">=15")) # first level will be treated as reference


## Age at last contact
CCSS_data$AGE_AT_LAST_CONTACT[CCSS_data$agelstcontact >= 0 & CCSS_data$agelstcontact < 25 ] <- "0-24"
CCSS_data$AGE_AT_LAST_CONTACT[CCSS_data$agelstcontact >= 25 & CCSS_data$agelstcontact < 35 ] <- "25-34"
CCSS_data$AGE_AT_LAST_CONTACT[CCSS_data$agelstcontact >= 35 & CCSS_data$agelstcontact < 45 ] <- "35-44"
CCSS_data$AGE_AT_LAST_CONTACT[CCSS_data$agelstcontact >= 45 ] <- ">=45"
CCSS_data$AGE_AT_LAST_CONTACT <- factor(CCSS_data$AGE_AT_LAST_CONTACT, levels = c("0-24", "25-34", "35-44", ">=45")) # first level will be treated as reference

## SN grade date
CCSS_data$gradedt <- as.numeric(CCSS_data$a_candx)

## datasets for subsequent neoplasm
subneo <- CCSS_data


# All unique samples in CCSS phenotype
PHENO.ANY_SN <- CCSS_data[c('ccssid', 'SEX', 'agedx', 'AGE_AT_DIAGNOSIS', 'agelstcontact', 
                           'chestrtgrp', 'neckrtgrp', 'abdomenrtgrp', 'abdomenrtgrp', 'brainrtgrp', 'pelvisrtgrp', 
                           'chestmaxrtdose', 'neckmaxrtdose', 'pelvismaxrtdose', 'abdmaxrtdose', 'maxsegrtdose', 'anth_DED5', 'alk_CED5', 'epipdose5', 'pt_cisED5')]
PHENO.ANY_SN <- PHENO.ANY_SN[!duplicated(PHENO.ANY_SN$ccssid),]

