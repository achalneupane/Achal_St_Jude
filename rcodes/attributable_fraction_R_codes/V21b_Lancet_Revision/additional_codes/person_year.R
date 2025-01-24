
## CCSS
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/")

ccss <- readRDS("CCSS_complete_data.rds")


# load("6.ccss_without_lifestyle.Any_SNs.Rdata")

length(unique(ccss$ccssid))
length(unique(PHENO.ANY_SN$ccssid))  

sum(PHENO.ANY_SN$PY)

############ CCSS #############
aa=table(ccss$ccssid)
length(aa)
### people with 2 ore more rows
more=aa[aa>1]
length(more)
### people with SN
ccss$age_sn=as.numeric(ccss$a_candx)
ccss_sn=ccss[ccss$a_candx>0,]
length(unique(ccss_sn$ccssid)) 
### if 2 or more rows in ccss_sn, then it is person with 2+ SNs.
aa=table(ccss_sn$ccssid)
length(aa)
more=aa[aa>1]
length(more)


ccss[ccss$ccssid==1000111,]


 library(dplyr)
ccss <- ccss %>% 
  group_by(ccssid) %>% 
  slice_head(n = 1)

### I assumed a_end is the age at the end of FU. For those with SN, a_candx is the age at SN. People were still at risk to develop 2nd or 3nrd.. SN after the 1st SN. So FU did not stop at a_candx but went to a_end.
ccss$age_end=as.numeric(ccss$a_end)
ccss$py=ccss$age_end-ccss$agedx-5
sum(ccss$py)  ### If "PHENO.ANY_SN" had the same # of people and same enty and end of FU, this should be the same as sum(PHENO.ANY_SN$PY). I see that they differ a bit. Achal should know the reason why, such as 7943 vs. 7918, and also for those with SN what is the end of FU (the last SN or till a_end)? Achal, I don't think there is a need to revise the analysis or change the results, either way they will be similar. 
# 184578.2 


## SJLIFE

sj <- readRDS("SJLIFE_complete_data.rds")

sj$gradedt <- as.Date(sj$gradedt, "%m/%d/%Y") ## **
sj$a_candx <- time_length(interval(as.Date(sj$dob), as.Date(sj$gradedt)), "years")

############ CCSS #############
aa=table(sj$sjlid)
length(aa)
### people with 2 ore more rows
more=aa[aa>1]
length(more)
### people with SN
sj$age_sn=as.numeric(sj$a_candx)
sj_sn=sj[sj$a_candx>0,]
length(unique(sj_sn$sjlid)) ##1636 with SN? paper said 1611
### if 2 or more rows in ccss_sn, then it is person with 2+ SNs.
aa=table(sj_sn$sjlid)
length(aa)
more=aa[aa>1]
length(more)




library(dplyr)
sj <- sj %>% 
  group_by(sjlid) %>% 
  slice_head(n = 1)

### I assumed a_end is the age at the end of FU. For those with SN, a_candx is the age at SN. People were still at risk to develop 2nd or 3nrd.. SN after the 1st SN. So FU did not stop at a_candx but went to a_end.
sj$age_end=as.numeric(sj$agelstcontact)
sj$py=sj$age_end-sj$agedx-5
sum(sj$py)  ### If "PHENO.ANY_SN" had the same # of people and same enty and end of FU, this should be the same as sum(PHENO.ANY_SN$PY). I see that they differ a bit. Achal should know the reason why, such as 7943 vs. 7918, and also for those with SN what is the end of FU (the last SN or till a_end)? Achal, I don't think there is a need to revise the analysis or change the results, either way they will be similar. 
# 90153.8 



################################
## Follow up regardless of SN ##
################################

## CCSS
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/")

ccss <- readRDS("CCSS_complete_data.rds")


# load("6.ccss_without_lifestyle.Any_SNs.Rdata")

length(unique(ccss$ccssid))

############ CCSS #############
aa=table(ccss$ccssid)
length(aa)
### people with 2 ore more rows
more=aa[aa>1]
length(more)
### people with SN
ccss$age_sn=as.numeric(ccss$a_candx)
ccss_sn=ccss[ccss$a_candx>0,]
length(unique(ccss_sn$ccssid)) 
### if 2 or more rows in ccss_sn, then it is person with 2+ SNs.
aa=table(ccss_sn$ccssid)
length(aa)
more=aa[aa>1]
length(more)


ccss[ccss$ccssid==1000111,]


ccss$age_end=as.numeric(ccss$a_end)

library(dplyr)
ccss <- ccss %>% 
  group_by(ccssid, age_end) %>% 
  slice(which.max(age_end))

### I assumed a_end is the age at the end of FU. For those with SN, a_candx is the age at SN. People were still at risk to develop 2nd or 3nrd.. SN after the 1st SN. So FU did not stop at a_candx but went to a_end.
ccss$py=ccss$age_end-ccss$agedx-5
sum(ccss$py)  ### If "PHENO.ANY_SN" had the same # of people and same enty and end of FU, this should be the same as sum(PHENO.ANY_SN$PY). I see that they differ a bit. Achal should know the reason why, such as 7943 vs. 7918, and also for those with SN what is the end of FU (the last SN or till a_end)? Achal, I don't think there is a need to revise the analysis or change the results, either way they will be similar. 
# 184578.2 


## SJLIFE

sj <- readRDS("SJLIFE_complete_data.rds")

# sj <- sj[sj$sjlid %in% PHENO.ANY_SN$sjlid,]

sj$gradedt <- as.Date(sj$gradedt, "%m/%d/%Y") ## **
sj$a_candx <- time_length(interval(as.Date(sj$dob), as.Date(sj$gradedt)), "years")

############ CCSS #############
aa=table(sj$sjlid)
length(aa)
### people with 2 ore more rows
more=aa[aa>1]
length(more)
### people with SN
sj$age_sn=as.numeric(sj$a_candx)
sj_sn=sj[sj$a_candx>0,]
length(unique(sj_sn$sjlid)) ##1636 with SN? paper said 1611
### if 2 or more rows in ccss_sn, then it is person with 2+ SNs.
aa=table(sj_sn$sjlid)
length(aa)
more=aa[aa>1]
length(more)



sj$age_end=as.numeric(sj$agelstcontact)
sj <- sj %>% 
  group_by(sjlid, age_end) %>% 
  slice(which.max(age_end))

### I assumed a_end is the age at the end of FU. For those with SN, a_candx is the age at SN. People were still at risk to develop 2nd or 3nrd.. SN after the 1st SN. So FU did not stop at a_candx but went to a_end.
sj$py=sj$age_end-sj$agedx-5
sum(sj$py)  ### If "PHENO.ANY_SN" had the same # of people and same enty and end of FU, this should be the same as sum(PHENO.ANY_SN$PY). I see that they differ a bit. Achal should know the reason why, such as 7943 vs. 7918, and also for those with SN what is the end of FU (the last SN or till a_end)? Achal, I don't think there is a need to revise the analysis or change the results, either way they will be similar. 
# 90153.8 


## Person year:
# SN          SJ          CCSS
# SNs         83990.13    172649
# SMNs        85488.43    179148
# Meningioma  88868.49    183024.8
# NMSCs       88513.76    179547.8
# Breast      89627.74    182791.3
# Thyroid     89362.05    183071.8
# Sarcoma     89857.59    182751.9

(605/83990.13)*1000
(463/85488.43)*1000
(149/88868.49)*1000
(251/88513.76)*1000
(76/89627.74)*1000
(87/89362.05)*1000
(33/89857.59)*1000

(1611/172649)*1000
(762/179148)*1000
(256/183024.8)*1000
(728/179547.8)*1000
(290/182791.3)*1000
(163/183071.8)*1000
(61/182751.9)*1000
