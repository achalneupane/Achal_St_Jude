setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2")


CCSS_data <- read.delim("ExportedCCSS_data.txt", header = T, sep = "\t")

BMI.PA.SMK.DRK <- read.delim("ExportedCCSS_BMI_PA_Smk_drink.txt", header = T, sep = "\t")

BMI.PA.SMK.DRK <- BMI.PA.SMK.DRK[c("ccssid", "a_base", "a_fu1", "a_fu2", "a_fu3", "a_fu2007", "a_fu5", "a_fu6",
                          "cbmi_0",  "cbmi_2", "cbmi_2007", "cbmi_5",
                          "t_eqmodfu2", "cdc_fu2", "t_eqmodfu5", "cdc_fu5", "t_eqmodfu6", "cdc_fu6", 
                          "smkcatb", "smkcatf2", "smkcatf07", "smkcatf5",
                          "riskyb", "riskyf7", "riskyfu5")]

colnames(BMI.PA.SMK.DRK) <- c("ccssid", "base.age", "fu1.age", "fu2.age", "fu3.age", "fu7.age", "fu5.age", "fu6.age",
                          "base.bmi",  "fu2.bmi", "fu7.bmi", "fu5.bmi",
                          "fu2.MET", "fu2.CDC", "fu5.MET", "fu5.CDC", "fu6.MET", "fu6.CDC", 
                          "base.smk", "fu2.smk", "fu7.smk", "fu5.smk",
                          "base.riskydrk", "fu7.riskydrk", "fu5.riskydrk")


# Make columns uniform for all follow-ups by adding missing columns
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


## Reshape BMI.PA.SMK.DRK to long format

# dw <- BMI.PA.SMK.DRK[1:2,]

cc <- reshape(BMI.PA.SMK.DRK, direction='long', 
        varying=c('base.age', 'base.bmi', "base.MET", 'base.CDC', "base.smk", "base.riskydrk",
                  'fu1.age', 'fu1.bmi', "fu1.MET", 'fu1.CDC', "fu1.smk", "fu1.riskydrk",
                  'fu2.age', 'fu2.bmi', "fu2.MET", 'fu2.CDC', "fu2.smk", "fu2.riskydrk",
                  'fu3.age', 'fu3.bmi', "fu3.MET", 'fu3.CDC', "fu3.smk", "fu3.riskydrk",
                  'fu7.age', 'fu7.bmi', "fu7.MET", 'fu7.CDC', "fu7.smk", "fu7.riskydrk",
                  'fu5.age', 'fu5.bmi', "fu5.MET", 'fu5.CDC', "fu5.smk", "fu5.riskydrk",
                  'fu6.age', 'fu6.bmi', "fu6.MET", 'fu6.CDC', "fu6.smk", "fu6.riskydrk"),
        timevar='var',
        times=c('base', 'fu1', 'fu2', 'fu3', 'fu7', 'fu5', 'fu6'),
        v.names=c('age', 'bmi', 'MET', 'CDC', 'smk', 'riskydrk'),
        idvar='ccssid')


cc[cc == "."] <- NA

cc$age <- as.numeric(cc$age)

## BMI
bmi_iid_dob_18 = subset(cc, age >= 18)
bmi_iid_dob_18 <- bmi_iid_dob_18[!is.na(bmi_iid_dob_18$bmi), ]
bmi_iid_dob_18_sorted = bmi_iid_dob_18[order(bmi_iid_dob_18$ccssid, bmi_iid_dob_18$age, decreasing = FALSE),]
bmi_iid_dob_18_uniq = bmi_iid_dob_18_sorted[!duplicated(bmi_iid_dob_18_sorted$ccssid),]

## MET
MET_iid_dob_18 = subset(cc, age >= 18)
MET_iid_dob_18 <- MET_iid_dob_18[!is.na(MET_iid_dob_18$MET), ]
MET_iid_dob_18_sorted = MET_iid_dob_18[order(MET_iid_dob_18$ccssid, MET_iid_dob_18$age, decreasing = FALSE),]
MET_iid_dob_18_uniq = MET_iid_dob_18_sorted[!duplicated(MET_iid_dob_18_sorted$ccssid),]

## smk
smk_iid_dob_18 = subset(cc, age >= 18)
smk_iid_dob_18 <- smk_iid_dob_18[!is.na(smk_iid_dob_18$smk), ]
smk_iid_dob_18_sorted = smk_iid_dob_18[order(smk_iid_dob_18$ccssid, smk_iid_dob_18$age, decreasing = FALSE),]
smk_iid_dob_18_uniq = smk_iid_dob_18_sorted[!duplicated(smk_iid_dob_18_sorted$ccssid),]

## drinking
drk_iid_dob_18 = subset(cc, age >= 18)
drk_iid_dob_18 <- drk_iid_dob_18[!is.na(drk_iid_dob_18$riskydrk), ]
drk_iid_dob_18_sorted = drk_iid_dob_18[order(drk_iid_dob_18$ccssid, drk_iid_dob_18$age, decreasing = FALSE),]
drk_iid_dob_18_uniq = drk_iid_dob_18_sorted[!duplicated(drk_iid_dob_18_sorted$ccssid),]


###################
## CCSS_original ##
###################

##############
## CCSS_exp ##
##############
ccss_exp.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/merged.dat.fam")
sum(CCSS_data$ccssid %in% ccss_exp.samples$V2)
# 3078
sum(ccss_exp.samples$V2 %in% CCSS_data$ccssid)
# 2936

CCSS_exp <- CCSS_data[CCSS_data$ccssid %in% ccss_exp.samples$V2,]



###########################################################
## Phenotype radiation and chemotherapy: create tertiles ##
###########################################################
## Create tertiles for ccss org and ccss exp separately
treatments <- CCSS_data[!duplicated(CCSS_data$ccssid),]


