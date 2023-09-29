rm(list=ls())

##############
## CCSS_exp ##
##############
ccss_exp.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/plink_data/merged.dat.fam")

##############
## CCSS_org ##
##############
ccss_org.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/plink_data/merged_batch2.fam")
ccss_org.samples$V1 <- sub("^(\\d+)_.*", "\\1", ccss_org.samples$V1)

ccss_samples <- c(ccss_exp.samples$V1, ccss_org.samples$V1)
overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/SJLIFE_WGS_samples_overlap_with_CCSS_org_SNP_samples.txt", header = T)

sum(ccss_samples %in% overlaps$ccssid)
ccss_samples <- ccss_samples[!(ccss_samples %in% overlaps$ccssid)]

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/CCSS_Data_from_Huiqi/RE__CCSS_phenotype_data2")

CCSS_data <- read.delim("ExportedCCSS_data.txt", header = T, sep = "\t", stringsAsFactors = F)
CCSS_data <- CCSS_data[CCSS_data$ccssid %in% ccss_samples,]

BMI.PA.SMK.DRK <- read.delim("ExportedCCSS_BMI_PA_Smk_drink.txt", header = T, sep = "\t", stringsAsFactors = F)

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

## Keep only those present in ccss WGS
BMI.PA.SMK.DRK <- BMI.PA.SMK.DRK[BMI.PA.SMK.DRK$ccssid %in% ccss_samples,]

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

BMI.PA.SMK.DRK[BMI.PA.SMK.DRK == "."] <- NA


# Create a vector of survey column names
surveys <- c("base", "fu1", "fu2", "fu3", "fu5", "fu6", "fu7")

# Initialize vectors to store merged data
merged_ccssid <- integer(0)
merged_age <- character(0)
merged_bmi <- character(0)
merged_smk <- character(0)
merged_riskydrk <- character(0)
merged_CDC <- character(0)

# Loop through each ID
for (i in 1:nrow(BMI.PA.SMK.DRK)) {
  id <- BMI.PA.SMK.DRK$ccssid[i]
  
  # Initialize variables to store data for the youngest adult survey
  youngest_age <- Inf
  youngest_bmi <- NA
  youngest_smk <- NA
  youngest_riskydrk <- NA
  youngest_CDC <- NA
  
  # Loop through each survey and find the youngest adult age
  for (survey in surveys) {
    age <- as.numeric(BMI.PA.SMK.DRK[i, paste0(survey, ".age")])
    if (!is.na(age) && age >= 18 && age < youngest_age) {
      youngest_age <- age
      youngest_bmi <- as.character(BMI.PA.SMK.DRK[i, paste0(survey, ".bmi")])
      youngest_smk <- as.character(BMI.PA.SMK.DRK[i, paste0(survey, ".smk")])
      youngest_riskydrk <- as.character(BMI.PA.SMK.DRK[i, paste0(survey, ".riskydrk")])
      youngest_CDC <- as.character(BMI.PA.SMK.DRK[i, paste0(survey, ".CDC")])
    }
  }
  
  # Append data to merged vectors
  merged_ccssid <- c(merged_ccssid, id)
  merged_age <- c(merged_age, youngest_age)
  merged_bmi <- c(merged_bmi, youngest_bmi)
  merged_smk <- c(merged_smk, youngest_smk)
  merged_riskydrk <- c(merged_riskydrk, youngest_riskydrk)
  merged_CDC <- c(merged_CDC, youngest_CDC)
}

# Create a new data frame with merged data
merged_data <- data.frame(
  ccssid = merged_ccssid,
  min_age = merged_age,
  min_bmi = merged_bmi,
  min_smk = merged_smk,
  min_riskydrk = merged_riskydrk,
  min_CDC = merged_CDC
)


dim(merged_data)
# 7943


count_same <- function(row) {
  sum(row == row[1], na.rm = TRUE)
}


test.1$same_count <- apply(test.1, 1, count_same)
table(test.1$same_count)

# Print the merged data frame
print(merged_data)
1000043