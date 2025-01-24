
setwd("R:\\Biostatistics\\Biostatistics2\\Qi\\QiCommon\\St Jude\\Yadav\\Achal")

load("6.sjlife_without_lifestyle.Any_SNs.Rdata")
 ls()
 dim(PHENO.ANY_SN)
 
 
 
######################################
## Attributable fraction of Any SNs ##
######################################
dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + 
                maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                EAS + AFR + 
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

##########################
## Get predicted values ##
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")


## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = T) # Overall
table(dat_all$event)



#############
## tx only ##
#############

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all

dat_tx$any_chemo_missing <- "No" # **
dat_tx$epitxn_dose_5.category = "None" ## **

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
## AF by tx
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
# af_by_tx <- round(af_by_tx,3)  # I would not do the rounding to get more accurate results for now and in 95% CI below
af_by_tx

## Male
N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Male"], na.rm = TRUE)
af_by_tx.male = (N_all.male - N_no_tx) / N_all.male
# af_by_tx.male <- round(af_by_tx.male,3)
af_by_tx.male

###################==== bootstraps starts here: random sample with replacement the ids from the original data, make new data based on the selected sample, repeat the above procedure 1000 times to the 1000 values of "af_by_tx"#########

library(dplyr)
lastrow=dat_all %>% 
  group_by(sjlid) %>% 
  slice(which.max(AGE_AT_LAST_CONTACT.cs1))
 ### Qi note: I think AGE_AT_LAST_CONTACT.cs1 is your age spline 1 and hence it is the originally age row. The above code take the last row in each person, so I can see how many unique ids.
all_ids=lastrow[,"sjlid"]


af_by_tx.save=NULL
af_by_tx.save=rbind(af_by_tx.save,af_by_tx)

af_by_tx.male.save=NULL
af_by_tx.male.save=rbind(af_by_tx.male.save,af_by_tx.male)


for(boot in 1:10){
#	boot=1;
print(boot)
set.seed(boot);
rows_sample=sample(1:nrow(all_ids), nrow(all_ids), replace = TRUE, prob = NULL)			  
length(rows_sample)
length(unique(rows_sample))
sample_ids=all_ids[rows_sample,]
##make fake ids -- consider each sampled id a new person
sample_ids$newid=seq(1:length(rows_sample))

#		 aa=table(sample_ids)
# table(aa)
### with boot=1 as an example, 1651 sjlids were seleted once, there was 4 ids were selected 6 times. For example SJL1264701
# aa[aa==6]

##### many to many merge to make new data
library(tidyverse)
df_test <- dat_all %>% 
  left_join(sample_ids, by = "sjlid")
  
#	df_test[1:1000,c("sjlid","start","end")]  
### with boot=1 as an example, SJL1264701 was selected 6 times, so there are 6 copies of this person's data in the new data.
df_test[df_test$sjlid=="SJL1264701",c("sjlid","start","end","newid")]


#### repeate the above to get AF with the bootstrapped data  ##########
fit_all = glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + 
                maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                EAS + AFR + 
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(df_test$PY), data = df_test)

df_test$pred_all = predict(fit_all, newdat = df_test, type = "response")
N_all = sum(df_test$pred_all, na.rm = T) # Overall

## Move relevant treatment exposures for everyone to no exposure
dat_tx = df_test

dat_tx$any_chemo_missing <- "No" # **
dat_tx$epitxn_dose_5.category = "None" ## **

df_test$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(df_test$pred_all, na.rm = TRUE)
N_no_tx = sum(df_test$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
af_by_tx



## Male
N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Male"], na.rm = TRUE)
af_by_tx.male = (N_all.male - N_no_tx) / N_all.male
# af_by_tx.male <- round(af_by_tx.male,3)
af_by_tx.male


af_by_tx.save=rbind(af_by_tx.save,af_by_tx)
af_by_tx.male.save=rbind(af_by_tx.male.save,af_by_tx.male)

}

af_by_tx.save[1,]  ## the original AF
#### use the 25th and 975th percentile as the 95% CI. #####
quantile(af_by_tx.save[-1,], probs = c(0.025,0.975))

##### In this exmaple with the overall AF, AF=0.0796835 with 95% CI 0.0533517 to 0.1060913 

af_by_tx = c(af_by_tx.save[1,], quantile(af_by_tx.save[-1,], probs = c(0.025,0.975))[1], quantile(af_by_tx.save[-1,], probs = c(0.025,0.975))[2])
af_by_tx.male = c(af_by_tx.male.save[1,], quantile(af_by_tx.male.save[-1,], probs = c(0.025,0.975))[1], quantile(af_by_tx.male.save[-1,], probs = c(0.025,0.975))[2])

