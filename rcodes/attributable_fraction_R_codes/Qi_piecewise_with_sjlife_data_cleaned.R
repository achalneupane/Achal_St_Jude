library(haven)
library(benchmarkme)
library(dplyr)
library(plyr)
library(data.table)
library (birk)
library(gtools)
library(stringr)
# library(tidyverse)
library(lubridate)
# benchmarkme::get_ram()
library(survival)

## function to add cubic spline
cubic_spline <- function(tp, knots)
{
  n=length(knots)
  f=rep(0,times=length(tp)*(n-1))
  dim(f)=c(length(tp),(n-1))
  #f1(x)=x#
  f[,1]=tp
  for(j in 1:(n-2))
  {
    for (i in 1:length(tp))
    {
      f[i,(j+1)]=max(0,(tp[i]-knots[j])^3, na.rm=TRUE)-max(0,(tp[i]-knots[n-1])^3, na.rm=TRUE)*(knots[n]-knots[j])/(knots[n]-knots[n-1])+
        max(0,(tp[i]-knots[n])^3, na.rm=TRUE)*(knots[n-1]-knots[j])/(knots[n]-knots[n-1])
    }
  }
  return(f)
}


data <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/toyadav.sas7bdat")
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
data <- cbind.data.frame(data, ethnicity.admixture[match(data$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


# data$agelstcontact.origina <- data$agelstcontact
data$gradedt <- data$evaldt ## SN grade date

## Age at SN diagnosis
data$AGE.ANY_SN <- time_length(interval(as.Date(data$dob), as.Date(data$gradedt)), "years")

## Attained age: age at last contact for cases is SN diagnosis date
data$agelstcontact[!is.na(data$evaldt)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]


## Add cubic spline of age attained 
breaks = seq(5, 95, 22.5)
cp = quantile(data$agelstcontact, breaks/100, na.rm = T)
cs = cubic_spline(data$agelstcontact, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
data <- cbind.data.frame(data, cs)


## define event (Any SN)
# 2= Meningioma, 4= Sarcoma, 5, Breast cancer, 6 = Thyroid, 7 = NMSC
data$event <- ifelse(data$sndx != 0, 1, 0) ## Any SN

### How many people had events?
length(unique(data$sjlid[data$event==1]))

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
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx and end is first event time; Achal: also added +5 in the line below
alldata$start[alldata$event==1 & alldata$first==1] <- alldata$agedx[alldata$event==1 & alldata$first==1] + 5
alldata$end[alldata$event==1 & alldata$first==1] <- as.numeric(difftime(alldata$gradedt[alldata$event==1 & alldata$first==1],alldata$dob[alldata$event==1 & alldata$first==1], units = 'days')/365.25)

#### For events that are not the first, segments are from the previous event date to this event date
alldata$previous <- as.Date(c(NA,alldata$gradedt[1:length(alldata$gradedt)-1]),origin = "1970-01-01")
gg1 <-cbind.data.frame(alldata$sjlid, alldata$event, alldata$first, alldata$previous, alldata$gradedt, alldata$start, alldata$end)
alldata[1:10,]
### if first>1 (i.e, 2 or more events, previous event time remained, others are missing)
alldata$previous[alldata$first==1] <- NA

alldata$start[alldata$first>1] <- as.numeric(difftime(alldata$previous[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)
alldata$end[alldata$first>1] <- as.numeric(difftime(alldata$gradedt[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)

### If one person has only 1 event, we need to have segment from agedx to event date (handled above), and then from event date to end of FU
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
# adata$end[adata$end<=adata$start] <- adata$end[adata$end<=adata$start] + 1/365
any <- adata[adata$end<=adata$start,] # 
diff=any$start-any$end 
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


SNs_py2=SNs_py
# SNs_py2=SNs_py[SNs_py$evt1==1,] 


SNs_py2$agedxcat <- factor(SNs_py2$agedxcat)
SNs_py2$sex <- factor(SNs_py2$sex)
SNs_py2$sex <- relevel(SNs_py2$sex, ref = "Male")
SNs_py2$braincat <- factor(SNs_py2$braincat)
SNs_py2$abdcat <- factor(SNs_py2$abdcat)
SNs_py2$chestcat <- factor(SNs_py2$chestcat)
SNs_py2$pelviscat <- factor(SNs_py2$pelviscat)
SNs_py2$neckcat <- factor(SNs_py2$neckcat)
SNs_py2$cat_anthra_jco_dose <- factor(SNs_py2$cat_anthra_jco_dose)
SNs_py2$cat_epitxn_dose <- factor(SNs_py2$cat_epitxn_dose)
SNs_py2$cat_aa_class_dose <- factor(SNs_py2$cat_aa_class_dose)



###############
## Model fit ##
###############
## 1. Any SN
fit_all <- glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                agedxcat + sex +
                braincat + abdcat + chestcat + cat_epitxn_dose +
                EAS + AFR,
               family = "poisson", offset = log(SNs_py2$PY), data = SNs_py2)

summary_fit_all <- summary(fit_all)


# Extract coefficients, standard errors, and p-values from the model summary
coefficients <- coef(summary_fit_all)
std_errors <- coefficients[, "Std. Error"]
p_values <- coefficients[, "Pr(>|z|)"]

# Calculate relative risks (exponentiated coefficients)
relative_risks <- exp(coefficients[, "Estimate"])

# Calculate 95% confidence intervals for relative risks
conf_int_low <- exp(coefficients[, "Estimate"] - 1.96 * std_errors)
conf_int_high <- exp(coefficients[, "Estimate"] + 1.96 * std_errors)

# Create a data frame to display the results
results.anySN <- data.frame(
  Coefficient = rownames(coefficients),
  Relative_Risk = relative_risks,
  CI_Lower = conf_int_low,
  CI_Upper = conf_int_high,
  P_Value = p_values
)

# Print the results
print(results.anySN)


