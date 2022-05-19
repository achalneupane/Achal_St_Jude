## Achal Neupane
## Date: 04/26/2022

################
## Event data ##
################

setwd("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data")
# install.packages("haven")
library(haven)
library(benchmarkme)
library(dplyr)
library(plyr)
library(data.table)
# benchmarkme::get_ram()
all.sas.files <- list.files()

##########################
## ctcaegrades.sas7bdat ##
##########################
data <- read_sas("ctcaegrades.sas7bdat")
data.df <- as.data.frame(data)
# View(data.df)

## Recode any values other than 0-5 as NAs
data.df$grade[!grepl("0|1|2|3|4|5", data.df$grade)] <- NA

## Omit NA grades
# data.df <- data.df[!is.na(data.df$grade),]

## Clean condition stings
data.df$condition <- gsub("_$", "", (gsub("_+", "_",  gsub("[^A-Za-z0-9]", "_", data.df$condition))))
dim(data.df)

# data.df <- data.df[c(2,5,6,8,9,10,11)]
data.df <- data.df[c(2,8,10,11)]

# > head(as.data.frame(data.df))
# sjlid            condition grade ageevent
# 1 SJL1527107         Hearing_Loss     1 58.83798
# 3 SJL1527107         Hearing_Loss     4 68.00000
# 4 SJL1527107         Hearing_Loss     4 70.00000
# 5 SJL1527107 Aortic_Root_Aneurysm     0 64.48456
# 6 SJL1527107 Aortic_Root_Aneurysm     0 70.96675
# 7 SJL1527107        Atrial_myxoma     0 64.48456

# data.df.condition <- data.df %>% group_by(condition) %>% top_n(1, ageevent)
SJLID <- unique(data.df$sjlid)



LIST2 <- list()
for(i in 1:length(SJLID)){
  tmp.df <- data.df[data.df$sjlid %in% SJLID[i],]
  CONDITIONS <- unique(tmp.df$condition)
  LIST1 <- list()
  for(j in 1:length(CONDITIONS)){
    print(paste0("Doing SJLID-- ", SJLID[i], " for condition-- ", CONDITIONS[j], " -- ITERATION i: ", i, " and j ", j))
    tmp.df.condition <- tmp.df[tmp.df$condition %in% CONDITIONS[j],]
    tmp.df.condition <- tmp.df.condition %>% arrange(condition, desc(grade), desc(ageevent)) %>% distinct(condition, .keep_all = TRUE)
    # tmp.df.condition <- as.data.frame(tmp.df.condition %>% group_by(condition) %>% top_n(1, grade) %>% top_n(1, ageevent))
    colnames(tmp.df.condition)[-c(1:2)] <- paste(tmp.df.condition$condition,colnames(tmp.df.condition)[-c(1:2)], sep = "_")
    LIST1[[j]] <- tmp.df.condition[3:4]
  }
  LIST2.tmp <- data.frame(LIST1)
  LIST2.tmp$sjlid <- SJLID[i]
  LIST2[[i]] <- LIST2.tmp
}

data.long <- rbindlist(LIST2, fill = TRUE)
# options(digits=15)
data.long <- data.long %>%
  relocate(sjlid)

# write.table(data.long, "Z:/ResearchHome/ClusterHome/aneupane/OUTPUT_FILES/ctcaegrades_Event_data_sjlife_horizontal_format_04_27_2022.csv", quote = FALSE, row.names = FALSE, sep = "\t")



# ###################################
# ## ctcaegrades_20210401.sas7bdat ##
# ###################################
# 
# 
# i = grep("ctcaegrades_20210401.sas7bdat", all.sas.files)
# data <- read_sas(all.sas.files[i])
# data.df <- as.data.frame(data)
# View(data.df)
# 
# ## Recode any values other than 0-5 as NAs
## Clean condition stings
# data.df$condition <- gsub("_$", "", (gsub("_+", "_",  gsub("[^A-Za-z0-9]", "_", data.df$condition))))
# dim(data.df)
# 
# ## Omit NA grades
# # data.df <- data.df[!is.na(data.df$grade),]
# 
# # data.df <- data.df[c(2,5,6,8,9,10,11)]
# data.df <- data.df[c(2,8,10,11)]
# 
# # > head(as.data.frame(data.df))
# # sjlid            condition grade ageevent
# # 1 SJL1527107         Hearing_Loss     1 58.83798
# # 3 SJL1527107         Hearing_Loss     4 68.00000
# # 4 SJL1527107         Hearing_Loss     4 70.00000
# # 5 SJL1527107 Aortic_Root_Aneurysm     0 64.48456
# # 6 SJL1527107 Aortic_Root_Aneurysm     0 70.96675
# # 7 SJL1527107        Atrial_myxoma     0 64.48456
# 
# # data.df.condition <- data.df %>% group_by(condition) %>% top_n(1, ageevent)
# SJLID <- unique(data.df$sjlid)
# 
# 
# 
# LIST2 <- list()
# for(i in 1:length(SJLID)){
#   tmp.df <- data.df[data.df$sjlid %in% SJLID[i],]
#   CONDITIONS <- unique(tmp.df$condition)
#   LIST1 <- list()
#   for(j in 1:length(CONDITIONS)){
#     print(paste0("Doing SJLID-- ", SJLID[i], " for condition-- ", CONDITIONS[j], " -- ITERATION i: ", i, " and j ", j))
#     tmp.df.condition <- tmp.df[tmp.df$condition %in% CONDITIONS[j],]
#     tmp.df.condition <-  tmp.df.condition <- tmp.df.condition %>% arrange(condition, desc(grade), desc(ageevent)) %>% distinct(condition, .keep_all = TRUE)
#     # tmp.df.condition <- as.data.frame(tmp.df.condition %>% group_by(condition) %>% top_n(1, grade) %>% top_n(1, ageevent))
#     colnames(tmp.df.condition)[-c(1:2)] <- paste(tmp.df.condition$condition,colnames(tmp.df.condition)[-c(1:2)], sep = "_")
#     LIST1[[j]] <- tmp.df.condition[3:4]
#   }
#   LIST2.tmp <- data.frame(LIST1)
#   LIST2.tmp$sjlid <- SJLID[i]
#   LIST2[[i]] <- LIST2.tmp
# }

# 
# data.long <- rbindlist(LIST2, fill = TRUE)
# # options(digits=10)
# data.long <- data.long %>%
#   relocate(sjlid)
# 
# write.table(data.long, "Z:/ResearchHome/ClusterHome/aneupane/OUTPUT_FILES/ctcaegrades_20210401_Event_data_sjlife_horizontal_format_Achal_04_26_2022.csv", quote = FALSE, row.names = FALSE, sep = "\t")


####################################################################################
####################################################################################
####################################################################################


##############
## Clinical ##
##############
setwd("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/")

##################
## Demographics ##
##################
data.demographics <- read_sas("demographics.sas7bdat")
data.demographics.df <- as.data.frame(data.demographics)
# View(data.demographics.df)
# cat(colnames(data.demographics.df), sep = ", ")

dim(data.demographics.df)
MRN <- data.demographics.df[1:2]
write.table(MRN, "Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/PHENOTYPE/MRN_SJLID_05_02_2022.csv", quote = FALSE, row.names = FALSE, sep = "\t")

data.demographics.df <- data.demographics.df[, c("sjlid", "gender", "race", "racegrp", "ethnic", "hispanic")]

###############
## Diagnosis ##
###############
data.diagnosis <- read_sas("diagnosis.sas7bdat")
data.diagnosis.df <- as.data.frame(data.diagnosis)
# View(data.diagnosis.df)
data.diagnosis.df <- data.diagnosis.df[, c("sjlid", "agedx", "diaggrp", "diag", "primdx")]

data.diagnosis.df <- data.diagnosis.df[data.diagnosis.df$primdx == 1 & !is.na(data.diagnosis.df$primdx),]

# data.diagnosis.df <- data.diagnosis.df %>% arrange(!is.na(diaggrp), desc(primdx)) %>%
#   distinct(sjlid, .keep_all = TRUE)
# 
# # excluding one sample with primdx == 0
# data.diagnosis.df <- data.diagnosis.df[data.diagnosis.df$primdx !=0,]
# # cat(colnames(data.diagnosis.df), sep = ", ")

data.diagnosis.df <- data.diagnosis.df[, c("sjlid", "agedx", "diaggrp", "diag")]

##################
## Last contact ##
##################
data.lastcondt <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Tracking Data/lstcondt.sas7bdat")
data.lastcondt.df <- as.data.frame(data.lastcondt)
View(data.lastcondt.df)

cat(colnames(data.lastcondt.df), sep = ", ")

data.lastcondt.df <- data.lastcondt.df[, c("sjlid", "agelstcontact")]

###################
## Chemosum_dose ##
###################
data.chemosum_dose <- read_sas("chemosum_dose.sas7bdat")
data.chemosum_dose.df <- as.data.frame(data.chemosum_dose)
View(data.chemosum_dose.df)

data.chemosum_dose.df <- data.chemosum_dose.df[, c("sjlid", "anthra_jama_dose_any", "anthra_jco_dose_any")]

###################
## Chemosum_yn ##
###################
data.chemosum_yn <- read_sas("chemosum_yn.sas7bdat")
data.chemosum_yn.df <- as.data.frame(data.chemosum_yn)
View(data.chemosum_yn.df)
cat(colnames(data.chemosum_yn.df), sep = ", ")

data.chemosum_yn.df <- data.chemosum_yn.df[, c("sjlid", "anthracyclines_any")]

#####################
## Radiationsum_yn ##
#####################
data.radiationsum_yn <- read_sas("radiationsum_yn.sas7bdat")
data.radiationsum_yn.df <- as.data.frame(data.radiationsum_yn)
View(data.radiationsum_yn.df)

# cat(colnames(data.radiationsum_yn.df), sep = ", ")

data.radiationsum_yn.df <- data.radiationsum_yn.df[colnames(data.radiationsum_yn.df)[grepl("^sjlid$|chest", colnames(data.radiationsum_yn.df), ignore.case = T)]]

#########################
## radiation_dosimetry ##
#########################
radiation_dosimetry <- read_sas("radiation_dosimetry.sas7bdat")
radiation_dosimetry.df <- as.data.frame(radiation_dosimetry)
View(radiation_dosimetry.df)

radiation_dosimetry.df <- radiation_dosimetry.df[, c("sjlid", "chestrt_yn", "maxchestrtdose")]

###################################
## rtdosimetrysjl_heart_20211112 ##
###################################
rtdosimetrysjl_heart_20211112 <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/MDA Dosimetry (Partial)/rtdosimetrysjl_heart_20211112.sas7bdat")
rtdosimetrysjl_heart_20211112.df <- as.data.frame(rtdosimetrysjl_heart_20211112)
View(rtdosimetrysjl_heart_20211112.df)

rtdosimetrysjl_heart_20211112.df <- rtdosimetrysjl_heart_20211112.df[, c("sjlid", "ccss", "RadYN", "HeartAvg")]



## Merge all dataframes
tt <- Reduce(
  function(x, y)
    merge(x, y, by = "sjlid", all = TRUE),
  list(data.demographics.df, data.diagnosis.df, data.lastcondt.df, 
       data.chemosum_dose.df, data.chemosum_yn.df, data.radiationsum_yn.df,
       radiation_dosimetry.df, rtdosimetrysjl_heart_20211112.df, data.long))

write.table(tt, "Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/PHENOTYPE/ctcaegrades_Event_data_sjlife_horizontal_format_05_02_2022.csv", quote = FALSE, row.names = FALSE, sep = "\t")


# count_dims <- function(x){
#   print(paste0("Unique: ", length(unique(x$sjlid)), " all: ", length(x$sjlid)))
# }
# 
# 
# mylist <- list(data.demographics.df, data.diagnosis.df, data.lastcondt.df, data.chemosum_dose.df, data.chemosum_yn.df, data.radiationsum_yn.df, radiation_dosimetry.df, rtdosimetrysjl_heart_20211112.df)
# names(mylist) <- c("data.demographics.df", "data.diagnosis.df", "data.lastcondt.df", "data.chemosum_dose.df", "data.chemosum_yn.df", "data.radiationsum_yn.df", "radiation_dosimetry.df", "rtdosimetrysjl_heart_20211112.df")
# lapply(mylist, count_dims)


