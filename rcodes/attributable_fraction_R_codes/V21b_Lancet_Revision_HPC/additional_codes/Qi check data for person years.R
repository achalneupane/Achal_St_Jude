

setwd("R:\\Biostatistics\\Biostatistics2\\Qi\\QiCommon\\St Jude\\Yadav\\Achal")

ccss <- readRDS("CCSS_complete_data.rds")
sj <- readRDS("SJLIFE_complete_data.rds")

load("6.ccss_without_lifestyle.Any_SNs.Rdata")

length(unique(ccss$ccssid)) ## 7943 unique people in demographic data
length(unique(PHENO.ANY_SN$ccssid))  ###7918 unique people in regression data  == what are the people removed? I will igore becuase the paper said N=7943 in CCSS.

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
length(unique(ccss_sn$ccssid)) ##1636 with SN? paper said 1611
### if 2 or more rows in ccss_sn, then it is person with 2+ SNs.
aa=table(ccss_sn$ccssid)
length(aa)
more=aa[aa>1]
length(more)


ccss[ccss$ccssid==1000111,]


 library(dplyr)
ccss=ccss %>% 
  group_by(ccssid) %>% 
  slice(which.max(ccssid))

### I assumed a_end is the age at the end of FU. For those with SN, a_candx is the age at SN. People were still at risk to develop 2nd or 3nrd.. SN after the 1st SN. So FU did not stop at a_candx but went to a_end.
ccss$age_end=as.numeric(ccss$a_end)
ccss$py=ccss$age_end-ccss$agedx-5
sum(ccss$py)  ### If "PHENO.ANY_SN" had the same # of people and same enty and end of FU, this should be the same as sum(PHENO.ANY_SN$PY). I see that they differ a bit. Achal should know the reason why, such as 7943 vs. 7918, and also for those with SN what is the end of FU (the last SN or till a_end)? Achal, I don't think there is a need to revise the analysis or change the results, either way they will be similar. 
 





