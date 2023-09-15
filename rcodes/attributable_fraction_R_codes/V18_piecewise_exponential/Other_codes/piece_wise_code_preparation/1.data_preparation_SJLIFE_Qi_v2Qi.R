library(lubridate)
library("survival")
library(haven)
library(dplyr)

rm(list=ls())
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")

dat1 <- PHENO.ANY_SN[c("sjlid", "dob", "agelstcontact", "agedx")]

## SN
subneo <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/subneo.sas7bdat")
deaths <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Tracking Data//deathsummary.sas7bdat")

head(subneo)
table(subneo$diaggrp)
dim(subneo)
# 1731 9

sum(subneo$sjlid %in% PHENO.ANY_SN$sjlid)
subneo <- subneo[subneo$sjlid %in% PHENO.ANY_SN$sjlid,] ## exclude those not in WGS

subneo$dob <- PHENO.ANY_SN$dob[match(subneo$sjlid, PHENO.ANY_SN$sjlid)]
# subneo$diagdt <- PHENO.ANY_SN$diagdt[match(subneo$sjlid, PHENO.ANY_SN$sjlid)]
subneo$agedx <- PHENO.ANY_SN$agedx[match(subneo$sjlid, PHENO.ANY_SN$sjlid)]

## Remove SN within 5 years of primary cancer
subneo$gradedt[grepl("3003-04-22", subneo$gradedt)] <- "2003-04-22"
subneo$AGE.ANY_SN <- time_length(interval(as.Date(subneo$dob), as.Date(subneo$gradedt)), "years")

## These two dates should be the (almost) same
# subneo$AGE.ANY_SN.after.childhood.cancer <- time_length(interval(as.Date(subneo$diagdt), as.Date(subneo$gradedt)), "years")
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx

# How many SNs after 5 years
subneo <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo$sjlid))
# 612


## For Any SN
subneo <- as.data.frame(subneo[c("sjlid", "gradedt", "diag", "diaggrp", "AGE.ANY_SN")])
# subneo$agedeath <- deaths$agedeath[match(subneo$sjlid, deaths$sjlid)]

# ## For breast cancer
# subneo <- subneo[grepl("breast", subneo$diag, ignore.case = T),]
merged_df <- merge(dat1, subneo, by = "sjlid", all = TRUE)
# View(merged_df)

# rm(list=ls()[!grepl("merged_df", ls())])
# ## De-identify IDs
# deidentify <- function(x) {
#   return(paste0("DE_Identified_", match(x, unique(x))))
# }

# merged_df$sjlid <- deidentify(merged_df$sjlid)

# # save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/merged_df.RData")
# 
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/merged_df.RData")

data <- merged_df[c("sjlid", "dob", "agelstcontact", "agedx", "gradedt", "AGE.ANY_SN")]

data <- data  %>% 
  distinct()  %>%
  arrange(sjlid, gradedt)

## Merge Phenotype to it
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_diet.Any_SNs.V14-4-3.Rdata")
data <- data[data$sjlid %in% PHENO.ANY_SN$sjlid,]

data <- cbind.data.frame(data, PHENO.ANY_SN[match(data$sjlid, PHENO.ANY_SN$sjlid), c("Pleiotropy_PRSWEB_PRS.tertile.category",
                                                                                     "AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4",
                                                                                     "AGE_AT_DIAGNOSIS", "gender", 
                                                                                     "maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category", "epitxn_dose_5.category",
                                                                                     "Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn",
                                                                                     "EAS", "AFR", 
                                                                                     "any_lifestyle_missing", "any_tx_missing")])




sum(data$AGE.ANY_SN > data$agelstcontact, na.rm = T) # e.g., SJL1391401
# SN.after.agelst <- cbind.data.frame(sjlid = data$sjlid[which(data$AGE.ANY_SN > data$agelstcontact)],
#                                     agelstcontact = data$agelstcontact[which(data$AGE.ANY_SN > data$agelstcontact)],
#                                     AGE.ANY_SN = data$AGE.ANY_SN[which(data$AGE.ANY_SN > data$agelstcontact)],
#                                     agediff = abs(data$agelstcontact[which(data$AGE.ANY_SN > data$agelstcontact)]- data$AGE.ANY_SN[which(data$AGE.ANY_SN > data$agelstcontact)]),
#                                     ageatdeath = data$agedeath[which(data$AGE.ANY_SN > data$agelstcontact)])



data$agelstcontact[which(data$AGE.ANY_SN > data$agelstcontact)] <- data$AGE.ANY_SN[which(data$AGE.ANY_SN > data$agelstcontact)]


data$event <- ifelse(!is.na(data$gradedt), 1, 0) # those with SN grade dates

rm(list = ls()[!ls() %in% "data"])

data$first <- ave(data$agelstcontact, data$sjlid, FUN = seq_along)
M <- max(data$first, na.rm = T)  ### maximum number of events
data[data$first==M,]$sjlid  #the id for the person with the maximum number of rows.
data[data$sjlid==data[data$first==M,]$sjlid,]

### Take the last row so we know the maximum number of events per person
event.number <- do.call(rbind, lapply(split(data, data$sjlid), tail, 1))[,c("sjlid","first")]
colnames(event.number) <- c("sjlid","maxE")

alldata <- merge(data,event.number,by.x="sjlid",by.y="sjlid")

alldata$start <- NULL
alldata$end <- NULL
### For those without event, start is agedx and end is Fu date === for SN, analysis starts from 5 years post DX (SNs within 5 years have been removed from the analysis)
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0] + 5
### Achal: Since the analysis start from 5 years post DX, the above line has been revised to: alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0]+5
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx and end is first event time
alldata$start[alldata$event==1 & alldata$first==1] <- alldata$agedx[alldata$event==1 & alldata$first==1] + 5
### Achal: also added +5 in the above line
alldata$end[alldata$event==1 & alldata$first==1] <- as.numeric(difftime(alldata$gradedt[alldata$event==1 & alldata$first==1],alldata$dob[alldata$event==1 & alldata$first==1], units = 'days')/365.25)

#### For events that are not the first, segments are from the previous event date to this event date
alldata$previous <- as.Date(c(NA,alldata$gradedt[1:length(alldata$gradedt)-1]),origin = "1970-01-01")
alldata[1:10,]
### if first>1 (i.e, 2 or more events, previous event time remained, others are missing)
alldata$previous[alldata$first==1] <- NA
alldata$start[alldata$first>1] <- as.numeric(difftime(alldata$previous[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)
alldata$end[alldata$first>1] <- as.numeric(difftime(alldata$gradedt[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)


### If one person has only 1 event, we need to have segment from agedx to event date (handled above), and then from event date to end of FU
### If one person has multiple events, we need to add segments from the last event to end of Fu.
row_add <- alldata[alldata$event==1 & alldata$first==alldata$maxE,] ## this includes (1) people with only 1 event (2) The last row for people with multiple events
row_add$start <- as.numeric(difftime(row_add$gradedt, row_add$dob, units = 'days')/365.25)
row_add$end <- row_add$agelstcontact;
#### remember to make these rows event as 0, as we are adding the time after the last event date.
row_add$event <- 0
### If we do first event analysis, we would not need these. I need to add an indicator here, in case I only need to get the first event analysis segments.
row_add$evt1 <- 0;

alldata$evt1 <- 0;
alldata$evt1[alldata$first==1] <- 1
adata <- rbind(alldata,row_add)
#### order by person and start date
adata <- adata[order(adata$sjlid, adata$start, decreasing = FALSE),]
table(adata$event)## Double check event numebr is correct

###any stop time <=start time
adata$end[adata$end<=adata$start] <- adata$end[adata$end<=adata$start] + 1/365
any <- adata[adata$end<=adata$start,] # 
### 2 people had the stroke date on the Fu date. In this case, either get rid of the 2 lines [i.e, no time is follow-up after the last event date], or add 1 day on end date of these 2 segments, assuming there were followed up 1 more day. Will not make much difference. 1 day out of 365 days is 0.0027.
diff=any$start-any$end ###Qi: These are people who had SN after the last contact date. Just wonder why this could happen. While it may not make the results differ, I wonder is there any reason to keep the events but change last contact date to be SN+1day? Depends on why there are SN after last contact date.
dim(adata)
final <- adata[adata$end>adata$start,]

minimum  <-  min(final$start, na.rm = TRUE) - 1
if (minimum < 0) minimum <- 0
maximum  <-  max(final$end, na.rm = TRUE) + 1
SNs_py <-  survSplit(final, cut=seq(minimum, maximum, 1), end="end",start="start",event="event") 
table(SNs_py$event)   ## Double check event numebr is correct
#### If you need the rows for first event analysis, take evt1=1
length(unique(SNs_py$sjlid[SNs_py$event==1]))
length(unique(SNs_py$sjlid[SNs_py$event==1 & SNs_py$evt1==1 ]))

SNs_py$PY <- SNs_py$end-SNs_py$start

###############
## Model fit ##
###############
SNs_py <- SNs_py[c("sjlid", "event", "Pleiotropy_PRSWEB_PRS.tertile.category",
  "AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4",
  "AGE_AT_DIAGNOSIS", "gender", 
  "maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category", "epitxn_dose_5.category", 
  "Current_smoker_yn", "PhysicalActivity_yn", "RiskyHeavyDrink_yn", "Obese_yn", 
  "EAS", "AFR", 
  "any_lifestyle_missing", "any_tx_missing", "PY","evt1")]

# SNs_py <- SNs_py[complete.cases(SNs_py),]

### Qi output the data to fit in SAS.
write.csv(SNs_py,file="R:/Biostatistics/Biostatistics2/Qi/QiCommon/St Jude/Achal/data.csv")

# library("geepack")
# fit_all <- geeglm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
#          AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#          AGE_AT_DIAGNOSIS + gender + 
#          maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category + 
#          Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn + 
#          EAS + AFR + 
#          any_lifestyle_missing + any_tx_missing,
#        family = "poisson", id = SNs_py$sjlid, offset = log(SNs_py$PY), corstr = "independence",  std.err = "san.se", data = SNs_py)
#####Qi: I tried the above model in SAS. If using corstr=cs in SAS, it had convergence issue. If I used type=ind then it converged in SAS. PRS not significant though, with p-value 0.17 and 0.14 below from SAS.
# Pleiotropy_PRSWEB_PR 2nd        0.2368   0.1729  -0.1021   0.5758    1.37   0.1709
# Pleiotropy_PRSWEB_PR 3rd        0.2401   0.1646  -0.0826   0.5628    1.46   0.1448
######## " AGE_AT_LAST_CONTACT.cs1", I  guess these are age cubic splines. One thing you can try is to make the age as categorical. For example, then one person would have many rows with the same age category. For example, one had 10 rows in 20-29 years old. Then you can collapse the rows into 1 row because the 10 rows have all X's the same. event=sum(event) and py=sum(py) over the 10 rows. In this way, you make your rows much smaller and will help with the convergence issue.


##### Qi: 	   
##### if not GEE (do first event analysis), the model works and here is how to get PAF.	   
SNs_py2=SNs_py[SNs_py$evt1==1,]
fit_all <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
			AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
         AGE_AT_DIAGNOSIS + gender + 
         maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category + 
         Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn + 
         EAS + AFR + 
         any_lifestyle_missing + any_tx_missing,
       family = "poisson", offset = log(SNs_py2$PY), data = SNs_py2)
#### to get PAF. 
#### First get the precited count from the original data and fit, the expected count should be the same as the # of events.	   
aa=predict(fit_all,type="response")   
sum(aa)	   

#### Assume you are interested in the PRS PAF. Make a fake new data, all people are in the reference of PRS.	   
pynew=SNs_py2
pynew$Pleiotropy_PRSWEB_PRS.tertile.category="1st"
table(pynew$Pleiotropy_PRSWEB_PRS.tertile.category)
bb=predict(fit_all,type="response",newdata=pynew)
PAF=(sum(aa)-sum(bb, na.rm = T))/sum(aa)  ##about 15%
PAF




