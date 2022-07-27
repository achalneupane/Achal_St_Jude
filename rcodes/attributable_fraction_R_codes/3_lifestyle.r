# #############################
# ## Adult habits/ Lifestyle ##
# #############################
# ## For each samples, get habits immediately after 18 years of age in agesurvey
# 
# # adlthabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.sas7bdat")
# adlthabits <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adlthabits.txt", header = T, sep ="\t")
# head(adlthabits)
# # remove duplicated rows
# adlthabits <- distinct(adlthabits)
# adlthabits$agesurvey <- floor(adlthabits$agesurvey_diff_of_dob_datecomp)
# 
# samples.sjlife <- unique(adlthabits$SJLIFEID)
# 
# lifestyle <- {}
# for (i in 1:length(samples.sjlife)){
#   print(paste0("Doing ", i))
#   dat <- adlthabits[adlthabits$SJLIFEID == samples.sjlife[i],]
#   if (max(dat$agesurvey) >= 18){
#     print("YES")
#     dat2 <- dat[dat$agesurvey >= 18,]
#     lifestyle.tmp <- dat2[which(dat2$agesurvey == min(dat2$agesurvey)),]
#     # } else {
#     #   print("NO")
#     #   lifestyle.tmp <-  dat[which(dat$agesurvey == max(dat$agesurvey)),]
#     #   lifestyle.tmp[9:ncol(lifestyle.tmp)] <- NA
#   }
#   lifestyle <- rbind.data.frame(lifestyle, lifestyle.tmp)
# }
# 
# sum(duplicated(lifestyle$SJLIFEID))
# lifestyle$SJLIFEID[duplicated(lifestyle$SJLIFEID)]
# ## Remove duplicate row
# lifestyle <- lifestyle[!duplicated(lifestyle$SJLIFEID),]
# 
# ## Add all samples
# # lifestyle <- cbind.data.frame(wgspop[,1:2], lifestyle[match(wgspop$MRN, lifestyle$mrn), ])
# # lifestyle <- lifestyle[-c(3,4)]
# # tt <- lifestyle
# ## Recode categorical variables
# lifestyle$relation[lifestyle$relation == 1] <- "Self"
# lifestyle$relation[lifestyle$relation == 2] <- "Parent"
# lifestyle$relation[lifestyle$relation == 3] <- "Other"
# 
# ## Recode smoker
# lifestyle$smoker[lifestyle$smoker == 1] <- "Past"
# lifestyle$smoker[lifestyle$smoker == 2] <- "Current"
# lifestyle$smoker[lifestyle$smoker == 3] <- "Never"
# lifestyle$smoker <- ifelse(lifestyle$smoker == "Current", 1, 0)
# 
# ## Recode to Y/N
# lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))][lifestyle[grepl("nopa|ltpa|bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))] == 1 ] <- 1
# lifestyle[grepl("nopa|ltpa", colnames(lifestyle))][lifestyle[grepl("nopa|ltpa", colnames(lifestyle))] == 2 ] <- 0
# lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))][lifestyle[grepl("bingedrink|heavydrink|heavydrink|riskydrink", colnames(lifestyle))] == 0 ] <- 0
# 
# # change the format of dates YYYY-MM-DD
# lifestyle$datecomp <- gsub("\\/", "-", lifestyle$datecomp)
# lifestyle$datecomp <- paste(sapply(strsplit(lifestyle$datecomp, "-"), `[`, 3), sapply(strsplit(lifestyle$datecomp, "-"), `[`, 1), sapply(strsplit(lifestyle$datecomp, "-"), `[`, 2), sep ="-")
# 
# lifestyle$dob <- gsub("\\/", "-", lifestyle$dob)
# lifestyle$dob <- paste(sapply(strsplit(lifestyle$dob, "-"), `[`, 3), sapply(strsplit(lifestyle$dob, "-"), `[`, 1), sapply(strsplit(lifestyle$dob, "-"), `[`, 2), sep ="-")
# 
# #######################
# ## Adolescent habits ##
# #######################
# 
# adolhabits <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adolhabits.sas7bdat")
# head(adolhabits)
# 
# ###############
# ## Adult BMI ##
# ###############
# 
# adultbmi <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/adultbmi.sas7bdat")
# head(adultbmi)
# 
# ## Add BMI and Nutrition to Lifestyle. Here, I am extracting the the BMI and Nutrition for the same age for each sample 
# lifestyle$BMI_KEY <- paste(lifestyle$SJLIFEID, lifestyle$agesurvey, sep = ":")
# 
# length(unique(adultbmi$sjlid))
# # 3640
# adultbmi$BMI_KEY <- paste(adultbmi$sjlid, adultbmi$AGE, sep = ":")
# sum(lifestyle$BMI_KEY %in% adultbmi$BMI_KEY)
# # 2964 
# ## samples that did not match by corresponding age
# cc <- lifestyle[!lifestyle$BMI_KEY %in% adultbmi$BMI_KEY,]
# 
# lifestyle <- cbind.data.frame(lifestyle, adultbmi[match(lifestyle$BMI_KEY, adultbmi$BMI_KEY), c("BMI", "HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")])
