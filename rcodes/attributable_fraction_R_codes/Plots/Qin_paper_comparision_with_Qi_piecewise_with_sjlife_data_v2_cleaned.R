rm(list=ls())


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
library(geepack)

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
data.all <- cbind.data.frame(data, ethnicity.admixture[match(data$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")

## read from datafreeze
chemo.2020 <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat')
chemo.2018 <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20181231/Clinical Data/chemosum_dose.sas7bdat')
# chemo.2020$match(PHENO.ANY_SN$sjlid, chemo.2020$sjlid)

radiation.2020 <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat')
radiation.2018 <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20181231/Clinical Data/radiation_dosimetry.sas7bdat')


data.all$maxsegrtdose.2020 <- radiation.2020$maxsegrtdose[match(data.all$sjlid, radiation.2020$sjlid)]
data.all$maxsegrtdose.2018 <- radiation.2018$maxsegrtdose[match(data.all$sjlid, radiation.2018$sjlid)]

data.all$maxchestrtdose.2020 <- radiation.2020$maxchestrtdose[match(data.all$sjlid, radiation.2020$sjlid)]
data.all$maxchestrtdose.2018 <- radiation.2018$maxchestrtdose[match(data.all$sjlid, radiation.2018$sjlid)]

data.all$maxneckrtdose.2020 <- radiation.2020$maxneckrtdose[match(data.all$sjlid, radiation.2020$sjlid)]
data.all$maxneckrtdose.2018 <- radiation.2018$maxneckrtdose[match(data.all$sjlid, radiation.2018$sjlid)]

data.all$maxpelvisrtdose.2020 <- radiation.2020$maxpelvisrtdose[match(data.all$sjlid, radiation.2020$sjlid)]
data.all$maxpelvisrtdose.2018 <- radiation.2018$maxpelvisrtdose[match(data.all$sjlid, radiation.2018$sjlid)]

data.all$maxabdrtdose.2020 <- radiation.2020$maxabdrtdose[match(data.all$sjlid, radiation.2020$sjlid)]
data.all$maxabdrtdose.2018 <- radiation.2018$maxabdrtdose[match(data.all$sjlid, radiation.2018$sjlid)]

data.all$aa_class_dose.2020 <- chemo.2020$alkylating_dose_5[match(data.all$sjlid, chemo.2020$sjlid)]
data.all$aa_class_dose.2018 <- chemo.2018$aa_class_dose_5[match(data.all$sjlid, chemo.2018$sjlid)]

data.all$anthra_jco_dose.2020 <- chemo.2020$anthracyclines_dose_5[match(data.all$sjlid, chemo.2020$sjlid)]
data.all$anthra_jco_dose.2018 <- chemo.2018$anthra_jco_dose_5[match(data.all$sjlid, chemo.2018$sjlid)]

data.all$epitxn_dose.2020 <- chemo.2020$epipodophyllotoxins_dose_5[match(data.all$sjlid, chemo.2020$sjlid)]
data.all$epitxn_dose.2018 <- chemo.2018$epitxn_dose_5[match(data.all$sjlid, chemo.2018$sjlid)]


## Check rt and chemo in Qi's data
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]

## MaxsegRT
PHENO.ANY_SN$maxsegrtdose.2020 <- radiation.2020$maxsegrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxsegrtdose.2018 <- radiation.2018$maxsegrtdose[match(PHENO.ANY_SN$sjlid, radiation.2018$sjlid)]
table(PHENO.ANY_SN$maxsegrtdose.2020 == PHENO.ANY_SN$maxsegrtdose) ## 2020
table(PHENO.ANY_SN$maxsegrtdose.2018 == PHENO.ANY_SN$maxsegrtdose) 

# In Qis data
PHENO.ANY_SN$maxsegrtdose.QI <- data.rtchemocheck$maxsegrtdose [match(PHENO.ANY_SN$sjlid, data.rtchemocheck$sjlid)]
table(PHENO.ANY_SN$maxsegrtdose.QI == PHENO.ANY_SN$maxsegrtdose.2018)


## Maxchest
PHENO.ANY_SN$maxchestrtdose.2020 <- radiation.2020$maxchestrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxchestrtdose.2018 <- radiation.2018$maxchestrtdose[match(PHENO.ANY_SN$sjlid, radiation.2018$sjlid)]
table(PHENO.ANY_SN$maxchestrtdose.2020 == PHENO.ANY_SN$maxchestrtdose) ## 2020
table(PHENO.ANY_SN$maxchestrtdose.2018 == PHENO.ANY_SN$maxchestrtdose) 

# In Qis data
PHENO.ANY_SN$maxchestrtdose.QI <- data.rtchemocheck$chestmaxrtdose[match(PHENO.ANY_SN$sjlid, data.rtchemocheck$sjlid)]
table(PHENO.ANY_SN$maxchestrtdose.QI == PHENO.ANY_SN$maxchestrtdose.2018)


## Maxneck
PHENO.ANY_SN$maxneckrtdose.2020 <- radiation.2020$maxneckrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxneckrtdose.2018 <- radiation.2018$maxneckrtdose[match(PHENO.ANY_SN$sjlid, radiation.2018$sjlid)]
table(PHENO.ANY_SN$maxneckrtdose.2020 == PHENO.ANY_SN$maxneckrtdose) ## 2020
table(PHENO.ANY_SN$maxneckrtdose.2018 == PHENO.ANY_SN$maxneckrtdose) 

# In Qis data
PHENO.ANY_SN$maxneckrtdose.QI <- data.rtchemocheck$neckmaxrtdose[match(PHENO.ANY_SN$sjlid, data.rtchemocheck$sjlid)]
table(PHENO.ANY_SN$maxneckrtdose.QI == PHENO.ANY_SN$maxneckrtdose.2018)


## Maxabdomen
PHENO.ANY_SN$maxabdrtdose.2020 <- radiation.2020$maxabdrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxabdrtdose.2018 <- radiation.2018$maxabdrtdose[match(PHENO.ANY_SN$sjlid, radiation.2018$sjlid)]
table(PHENO.ANY_SN$maxabdrtdose.2020 == PHENO.ANY_SN$maxabdrtdose) ## 2020
table(PHENO.ANY_SN$maxabdrtdose.2018 == PHENO.ANY_SN$maxabdrtdose) 
PHENO.ANY_SN$MRN[which(PHENO.ANY_SN$maxabdrtdose.2020 != PHENO.ANY_SN$maxabdrtdose)]

## Maxpelvis
PHENO.ANY_SN$maxpelvisrtdose.2020 <- radiation.2020$maxpelvisrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxpelvisrtdose.2018 <- radiation.2018$maxpelvisrtdose[match(PHENO.ANY_SN$sjlid, radiation.2018$sjlid)]
table(PHENO.ANY_SN$maxpelvisrtdose.2020 == PHENO.ANY_SN$maxpelvisrtdose) ## 2020
table(PHENO.ANY_SN$maxpelvisrtdose.2018 == PHENO.ANY_SN$maxpelvisrtdose) 
PHENO.ANY_SN$MRN[which(PHENO.ANY_SN$maxpelvisrtdose.2020 != PHENO.ANY_SN$maxpelvisrtdose)]


## aa_class_dose_5 **
PHENO.ANY_SN$aa_class_dose_5.2020 <- chemo.2020$alkylating_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2020$sjlid)]
PHENO.ANY_SN$aa_class_dose_5.2018 <- chemo.2018$aa_class_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2018$sjlid)]
table(PHENO.ANY_SN$aa_class_dose_5.2020 == PHENO.ANY_SN$aa_class_dose_5) ## 2020
table(PHENO.ANY_SN$aa_class_dose_5.2018 == PHENO.ANY_SN$aa_class_dose_5) 
PHENO.ANY_SN$MRN[which(PHENO.ANY_SN$aa_class_dose_5.2020 != PHENO.ANY_SN$aa_class_dose_5)]

## anthra_jco_dose_5 **
PHENO.ANY_SN$anthra_jco_dose_5.2020 <- chemo.2020$anthracyclines_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2020$sjlid)]
PHENO.ANY_SN$anthra_jco_dose_5.2018 <- chemo.2018$anthra_jco_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2018$sjlid)]
table(PHENO.ANY_SN$anthra_jco_dose_5.2020 == PHENO.ANY_SN$anthra_jco_dose_5) ## 2020
table(PHENO.ANY_SN$anthra_jco_dose_5.2018 == PHENO.ANY_SN$anthra_jco_dose_5) 
PHENO.ANY_SN$MRN[which(PHENO.ANY_SN$anthra_jco_dose_5.2020 != PHENO.ANY_SN$anthra_jco_dose_5)]

## epitxn_dose_5 **
PHENO.ANY_SN$epitxn_dose_5.2020 <- chemo.2020$epipodophyllotoxins_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2020$sjlid)]
PHENO.ANY_SN$epitxn_dose_5.2018 <- chemo.2018$epitxn_dose_5[match(PHENO.ANY_SN$sjlid, chemo.2018$sjlid)]
table(PHENO.ANY_SN$epitxn_dose_5.2020 == PHENO.ANY_SN$epitxn_dose_5) ## 2020
table(PHENO.ANY_SN$epitxn_dose_5.2018 == PHENO.ANY_SN$epitxn_dose_5) 
PHENO.ANY_SN$MRN[which(PHENO.ANY_SN$epitxn_dose_5.2020 != PHENO.ANY_SN$epitxn_dose_5)]


## 1.-----------------------------------maxsegRT
table(data.all$braincat)
# 0    1    2    3    4 
# 3452   25  978  368  351

duplicated(data.all$sjlid)


data.all$braincat.AN <- cut(data.all$maxsegrtdose, breaks = c(0, 200, 1799, 2999, max(data.all$maxsegrtdose, na.rm = T)),
    labels = c(0:3),
    include.lowest = TRUE)
levels(data.all$braincat.AN) <- c(levels(data.all$braincat.AN), 4)
data.all$braincat.AN[is.na(data.all$braincat.AN)] <- 4

data.rtchemocheck$braincat.AN <- cut(data.rtchemocheck$maxsegrtdose, breaks = c(0, 200, 1799, 2999, max(data.rtchemocheck$maxsegrtdose, na.rm = T)),
                            labels = c(0:3),
                            include.lowest = TRUE)
levels(data.rtchemocheck$braincat.AN) <- c(levels(data.rtchemocheck$braincat.AN), 4)
data.rtchemocheck$braincat.AN[is.na(data.rtchemocheck$braincat.AN)] <- 4


data.rtchemocheck$braincat.AN.2018 <- cut(data.rtchemocheck$maxsegrtdose.2018, breaks = c(0, 200, 1799, 2999, max(data.rtchemocheck$maxsegrtdose.2018, na.rm = T)),
                                     labels = c(0:3),
                                     include.lowest = TRUE)
levels(data.rtchemocheck$braincat.AN.2018) <- c(levels(data.rtchemocheck$braincat.AN.2018), 4)
data.rtchemocheck$braincat.AN.2018[is.na(data.rtchemocheck$braincat.AN.2018)] <- 4


table(data.all$braincat.AN)
# 0    1    2    3    4 
# 3222   19  978  368  587 

table(data.rtchemocheck$braincat.AN)
table(data.rtchemocheck$braincat.AN.2018)
table(PHENO.ANY_SN$maxsegrtdose.category)

table(data.rtchemocheck$maxabdrtdose.2018 == data.rtchemocheck$maxsegrtdose)

cc <- data.all[c("sjlid", "braincat", "braincat.AN", "maxsegrtdose", "maxabdrtdose.2018", "maxabdrtdose.2020", "braindose")]
## There are NAs but not labelled 4 in Qi's data
# 0    1    2    3    4 
# 3222   19  978  368  587
# Also, for example, IDs SJL4186407 and SJL5056617 are incorrectly categorized; and when there are missing values in SJL5244907 or SJL5245607..


## Our data
table(PHENO.ANY_SN$maxsegrtdose.category)
# None   >0-<18 >=18-<30     >=30  Unknown 
# 3255       18      568      508       52 

table(data.rtchemocheck$braincat)
# 0    1    2    3    4 
# 3218   18  556  279  331 

table(data.rtchemocheck$braincat.AN)



## 2.-----------------------------------maxneckRT
table(data.all$neckcat)
# 0    1    2    3 
# 3801  414  530  429 

data.all$neckcat.AN <- cut(data.all$neckmaxrtdose, breaks = c(0, 200, 1099, 1999, 2999, max(data.all$neckmaxrtdose, na.rm = T)),
    labels = c(0:4),
    include.lowest = TRUE)
levels(data.all$neckcat.AN) <- c(levels(data.all$neckcat.AN), 5)
data.all$neckcat.AN[is.na(data.all$neckcat.AN)] <- 5
table(data.all$neckcat.AN)
# 0    1    2    3    4    5 
# 3652    5   77  504  354  582 


data.rtchemocheck$neckcat.AN.2020 <- cut(data.rtchemocheck$neckmaxrtdose, breaks = c(0, 200, 1099, 1999, 2999, max(data.rtchemocheck$neckmaxrtdose, na.rm = T)),
                           labels = c(0:4),
                           include.lowest = TRUE)
levels(data.rtchemocheck$neckcat.AN.2020) <- c(levels(data.rtchemocheck$neckcat.AN.2020), 5)
data.rtchemocheck$neckcat.AN.2020[is.na(data.rtchemocheck$neckcat.AN.2020)] <- 5
table(data.rtchemocheck$neckcat.AN.2020)

data.rtchemocheck$neckcat.AN.2018 <- cut(data.rtchemocheck$neckmaxrtdose, breaks = c(0, 200, 1099, 1999, 2999, max(data.rtchemocheck$neckmaxrtdose, na.rm = T)),
                                         labels = c(0:4),
                                         include.lowest = TRUE)
levels(data.rtchemocheck$neckcat.AN.2018) <- c(levels(data.rtchemocheck$neckcat.AN.2018), 5)
data.rtchemocheck$neckcat.AN.2018[is.na(data.rtchemocheck$neckcat.AN.2018)] <- 5
table(data.rtchemocheck$neckcat.AN.2018)


table(PHENO.ANY_SN$maxneckrtdose.category)
table(data.rtchemocheck$neckcat)
table(data.rtchemocheck$neckcat.AN.2020)
table(data.rtchemocheck$neckcat.AN.2018)

cc <- data.all[c("sjlid", "neckcat", "neckcat.AN", "neckmaxrtdose", "maxneckrtdose.2018", "maxneckrtdose.2020")]

## Our data
table(PHENO.ANY_SN$maxneckrtdose.category)
# None   >0-<11 >=11-<20 >=20-<30     >=30  Unknown 
# 3561        6       90      426      270       48 

## 3.-----------------------------------abdomenRT
table(data.all$abdcat)
# 0    1    2    3 
# 4039  567  309  259 

data.all$abdcat.AN <- cut(data.all$abdmaxrtdose, breaks = c(0, 200, 2999, max(data.all$abdmaxrtdose, na.rm = T)),
    labels = c(0:2),
    include.lowest = TRUE)
levels(data.all$abdcat.AN) <- c(levels(data.all$abdcat.AN), 3)
data.all$abdcat.AN [is.na(data.all$abdcat.AN)] <- 3
table(data.all$abdcat.AN)
# 0    1    2    3 
# 3718  563  310  583

## Our data
table(PHENO.ANY_SN$maxabdrtdose.category)
# None  >0-<30    >=30 Unknown 
# 3568     522     263      48


## 4.-----------------------------------chestRT
table(data.all$chestcat)
# 0    1    2    3 
# 3872  124  851  327 
data.all$chestcat.AN <- cut(data.all$chestmaxrtdose, breaks = c(0, 200, 1999, max(data.all$chestmaxrtdose, na.rm = T)),
                                            labels = c(0:2),
                                            include.lowest = TRUE)
levels(data.all$chestcat.AN) <- c(levels(data.all$chestcat.AN), 3)
data.all$chestcat.AN [is.na(data.all$chestcat.AN)] <- 3
table(data.all$chestcat.AN)
# 0    1    2    3 
# 3621  119  851  583

## Our data
table(PHENO.ANY_SN$maxchestrtdose.category)
# None  >0-<20    >=20 Unknown 
# 3498     149     705      49 


## 5.-----------------------------------pelvisRT
table(data.all$pelviscat)
# 0    1    2    3 
# 4206  145  597  226 

data.all$pelviscat.AN <- cut(data.all$pelvismaxrtdose, breaks = c(0, 200, 1999, max(data.all$pelvismaxrtdose, na.rm = T)),
                                             labels = c(0:2),
                                             include.lowest = TRUE)
levels(data.all$pelviscat.AN) <- c(levels(data.all$pelviscat.AN), 3)
data.all$pelviscat.AN [is.na(data.all$pelviscat.AN)] <- 3
table(data.all$pelviscat.AN)
# 0    1    2    3 
# 3853  141  597  583 
## Our data
table(PHENO.ANY_SN$maxpelvisrtdose.category)
# None  >0-<20    >=20 Unknown 
# 3669     169     518      45




## 6.-----------------------------------aa_class_dose
table(data.all$cat_aa_class_dose)
# 0    1    2    3 
# 2114  924  936 1159 

TERT = unname(quantile(data.rtchemocheck$aa_class_dose[data.rtchemocheck$aa_class_dose !=0], c(1/3, 2/3, 1), na.rm = T))
data.all$cat_aa_class_dose.AN <- cut(data.all$aa_class_dose, breaks = c(-Inf, 0, TERT),
                                               labels = c("None", "1st", "2nd", "3rd"),
                                               include.lowest = TRUE)

data.all$cat_aa_class_dose.AN <- cut(data.all$aa_class_dose., breaks = c(-Inf, 0, TERT),
                                     labels = c("None", "1st", "2nd", "3rd"),
                                     include.lowest = TRUE)

table(data.all$cat_aa_class_dose.AN)
# None  1st  2nd  3rd 
# 2114  926  933 1160 

table(PHENO.ANY_SN$aa_class_dose_5.category)
# None  1st  2nd  3rd 
# 1943  819  812  815
table(data.rtchemocheck$cat_aa_class_dose)
# 0    1    2    3 
# 1911  820  824  820 




data.all$aa_class_dose.AN <- PHENO.ANY_SN$aa_class_dose_5.2020[match(data.all$sjlid, PHENO.ANY_SN$sjlid)]
data.all$aa_class_dose_5.category <- PHENO.ANY_SN$aa_class_dose_5.category[match(data.all$sjlid, PHENO.ANY_SN$sjlid)]
table(data.all$aa_class_dose.AN == data.all$aa_class_dose)

table(PHENO.ANY_SN$aa_class_dose_5.category)
# None  1st  2nd  3rd 
# 2146  923  923 1155

table(data.rtchemocheck$cat_aa_class_dose)
table(PHENO.ANY_SN$aa_class_dose_5.category)

cc <- data.all[c("sjlid", "aa_class_dose", "aa_class_dose.AN", "cat_aa_class_dose", "cat_aa_class_dose.AN", "aa_class_dose_5.category")]

## 7.-----------------------------------anthra_jco_dose
table(data.all$cat_anthra_jco_dose)
# 0    1    2    3 
# 2379  965  911  894 

TERT = unname(quantile(data.rtchemocheck$anthra_jco_dose[data.rtchemocheck$anthra_jco_dose !=0], c(1/3, 2/3, 1), na.rm = T))
data.all$cat_anthra_jco_dose.AN <- cut(data.all$anthra_jco_dose, breaks = c(-Inf, 0, TERT),
                                     labels = c("None", "1st", "2nd", "3rd"),
                                     include.lowest = TRUE)
table(data.all$cat_anthra_jco_dose.AN )
# None  1st  2nd  3rd 
# 2379  924  932  914 

data.all$anthra_jco_dose.AN <- PHENO.ANY_SN$anthra_jco_dose_5[match(data.all$sjlid, PHENO.ANY_SN$sjlid)]
table(data.all$anthra_jco_dose.AN == data.all$anthra_jco_dose)
TERT = unname(quantile(data.rtchemocheck$anthra_jco_dose.AN[data.rtchemocheck$anthra_jco_dose.AN !=0], c(1/3, 2/3, 1), na.rm = T))
data.all$cat_anthra_jco_dose.AN <- cut(data.all$anthra_jco_dose., breaks = c(-Inf, 0, TERT),
                                       labels = c("None", "1st", "2nd", "3rd"),
                                       include.lowest = TRUE)

table(data.all$cat_anthra_jco_dose.AN.2020)
# None  1st  2nd  3rd 
# 2413  917  916  917

cc <- data.all[c("sjlid", "anthra_jco_dose", "anthra_jco_dose.AN", "cat_anthra_jco_dose", "cat_anthra_jco_dose.AN", "cat_anthra_jco_dose.AN.2020")]

table(data.rtchemocheck$cat_anthra_jco_dose)
table(PHENO.ANY_SN$anthra_jco_dose_5.category)



## 7.-----------------------------------epitxn_dose
table(data.all$cat_epitxn_dose)
# 0    1    2    3 
# 3521  574  529  528 
TERT = unname(quantile(data.rtchemocheck$epitxn_dose[data.rtchemocheck$epitxn_dose !=0], c(1/3, 2/3, 1), na.rm = T))
data.all$cat_epitxn_dose.AN <- cut(data.all$epitxn_dose, breaks = c(-Inf, 0, TERT),
                                       labels = c("None", "1st", "2nd", "3rd"),
                                       include.lowest = TRUE)
table(data.all$cat_epitxn_dose.AN)
# None  1st  2nd  3rd 
# 3521  544  543  544

data.all$epitxn_dose.AN <- PHENO.ANY_SN$epitxn_dose_5.2020[match(data.all$sjlid, PHENO.ANY_SN$sjlid)]
table(data.all$aa_class_dose.AN == data.all$aa_class_dose)


data.all$agedx
data.all$agedxcat.AN <- data.all$agedx
data.all$agedxcat.AN[data.all$agedx >= 0 & data.all$agedx < 5 ] <- "0-4"
data.all$agedxcat.AN[data.all$agedx >= 5 & data.all$agedx < 10 ] <- "5-9"
data.all$agedxcat.AN[data.all$agedx >= 10 & data.all$agedx < 15 ] <- "10-14"
data.all$agedxcat.AN[data.all$agedx >= 15 ] <- ">=15"
data.all$agedxcat.AN <- factor(data.all$agedxcat.AN, levels = c("0-4", "5-9", "10-14", ">=15")) # first level will be treated as reference
table(data.all$agedxcat.AN)
table(data.all$agedxcat)


############
## Any SN ##
############

data <- data.all
# data$agelstcontact.original <- data$agelstcontact
data$gradedt <- data$evaldt ## SN grade date

## Age at SN diagnosis
data$AGE.ANY_SN <- time_length(interval(as.Date(data$dob), as.Date(data$gradedt)), "years")

## Attained age: age at last contact for cases is SN diagnosis date
data$agelstcontact[!is.na(data$evaldt)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

## define event (Any SN)
# 2= Meningioma, 4= Sarcoma, 5, Breast cancer, 6 = Thyroid, 7 = NMSC
data$event <- ifelse(data$sndx != 0, 1, 0) ## Any SN

## This is when we limit to first event
## Check how many within 5 years of primary cancer
data$AGE.ANY_SN.after.childhood.cancer.from.agedx <- data$AGE.ANY_SN - data$agedx
sum(data$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5, na.rm = T)
# 0 # Since it's 0, they all are after 5 years of primary diagnosis
# Get first event after 5 years of primary diagnosis
table(data$event == 1)
CA <- data[data$sndx != 0,]
CO <- data[data$sndx == 0,]
CA <- setDT(CA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

data <- rbind.data.frame(CA, CO)


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

alldata$start <- NA
alldata$end <- NA
### For those without event, start is agedx and end is Fu date === for SN, analysis starts from 5 years post DX (SNs within 5 years have been removed from the analysis)
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0] + 5
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx +5 and end is first event time; Achal: also added +5 in the line below
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


# SNs_py2=SNs_py
SNs_py2=SNs_py[SNs_py$evt1==1,]


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



## Add cubic spline of age attained
breaks = seq(5, 95, 22.5)
cp = quantile(SNs_py2$end, breaks/100, na.rm = T)
cs = cubic_spline(SNs_py2$end, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
SNs_py2 <- cbind.data.frame(SNs_py2, cs)

# library(splines)
# # Define evenly spaced knots
# num_knots <- 5
# knots <- seq(min(SNs_py2$end), max(SNs_py2$end), length.out = num_knots + 2)[2:(num_knots + 1)]
# 
# # Create a cubic spline with 5 evenly spaced knots
# cubic_spline <- ns(SNs_py2$end, knots = knots)
# 
# colnames(cubic_spline) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4", "AGE_AT_LAST_CONTACT.cs5", "AGE_AT_LAST_CONTACT.cs6")
# cs <- cubic_spline
# SNs_py2 <- cbind.data.frame(SNs_py2, cs)
###############
## Model fit ##
###############
## 1. Any SN
fit_all <- glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                agedxcat + sex +
                braincat.AN + abdcat.AN + chestcat.AN + cat_epitxn_dose.AN +
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
  Relative_Risk = round(relative_risks,1),
  CI_Lower = round(conf_int_low,1),
  CI_Upper = round(conf_int_high,1),
  P_Value = p_values
)

# Print the results
print(results.anySN)

################
## Meningioma ##
################
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

data <- data.all

# data$agelstcontact.original <- data$agelstcontact
data$gradedt <- data$evaldt ## SN grade date

## Age at SN diagnosis
data$AGE.ANY_SN <- time_length(interval(as.Date(data$dob), as.Date(data$gradedt)), "years")

## Attained age: age at last contact for cases is SN diagnosis date
data$agelstcontact[!is.na(data$evaldt)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

## define event (Any SN)
# 2= Meningioma, 4= Sarcoma, 5, Breast cancer, 6 = Thyroid, 7 = NMSC
# data$event <- ifelse(data$sndx != 0, 1, 0) ## Any SN
data$event <- ifelse(data$sndx == 2, 1, 0) ## Meningioma

## This is when we limit to first event
## Check how many within 5 years of primary cancer
data$AGE.ANY_SN.after.childhood.cancer.from.agedx <- data$AGE.ANY_SN - data$agedx
sum(data$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5, na.rm = T)
# 0 # Since it's 0, they all are after 5 years of primary diagnosis
## Get first event after 5 years of primary diagnosis
table(data$event == 1)
CA <- data[data$sndx == 2,]
CO <- data[data$sndx == 0,]
# CA <- setDT(CA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

data <- rbind.data.frame(CA, CO)
dim(data)

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

alldata$start <- NA
alldata$end <- NA
### For those without event, start is agedx and end is Fu date === for SN, analysis starts from 5 years post DX (SNs within 5 years have been removed from the analysis)
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0] + 5
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx +5 and end is first event time; Achal: also added +5 in the line below
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



## Add cubic spline of age attained
breaks = seq(5, 95, 22.5)
cp = quantile(SNs_py2$end, breaks/100, na.rm = T)
cs = cubic_spline(SNs_py2$end, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
SNs_py2 <- cbind.data.frame(SNs_py2, cs)


# # library(splines)
# # # Define evenly spaced knots
# num_knots <- 4
# knots <- seq(min(SNs_py2$end), max(SNs_py2$end), length.out = num_knots + 2)[2:(num_knots + 1)]
# 
# # Create a cubic spline with 5 evenly spaced knots
# cubic_spline <- ns(SNs_py2$end, knots = knots)
# 
# colnames(cubic_spline) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4", "AGE_AT_LAST_CONTACT.cs5")
# # colnames(cubic_spline) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4", "AGE_AT_LAST_CONTACT.cs5", "AGE_AT_LAST_CONTACT.cs6")
# cs <- cubic_spline
# SNs_py2 <- cbind.data.frame(SNs_py2, cs)

###############
## Model fit ##
###############
## 2. Meningioma
fit_all <- glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                 agedxcat + sex +
                 braincat + cat_epitxn_dose,
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
results.Meningioma <- data.frame(
  Coefficient = rownames(coefficients),
  Relative_Risk = round(relative_risks,1),
  CI_Lower = round(conf_int_low,1),
  CI_Upper = round(conf_int_high,1),
  P_Value = p_values
)

# Print the results
print(results.Meningioma)

#############
## Sarcoma ##
#############
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

data <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/toyadav.sas7bdat")
ethnicity.admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
data <- cbind.data.frame(data, ethnicity.admixture[match(data$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


# data$agelstcontact.original <- data$agelstcontact
data$gradedt <- data$evaldt ## SN grade date

## Age at SN diagnosis
data$AGE.ANY_SN <- time_length(interval(as.Date(data$dob), as.Date(data$gradedt)), "years")

## Attained age: age at last contact for cases is SN diagnosis date
data$agelstcontact[!is.na(data$evaldt)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

## define event (Any SN)
# 2= Meningioma, 4= Sarcoma, 5, Breast cancer, 6 = Thyroid, 7 = NMSC
data$event <- ifelse(data$sndx == 4, 1, 0) ## Sarcoma

## This is when we limit to first event
## Check how many within 5 years of primary cancer
data$AGE.ANY_SN.after.childhood.cancer.from.agedx <- data$AGE.ANY_SN - data$agedx
sum(data$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5, na.rm = T)
# 0 # Since it's 0, they all are after 5 years of primary diagnosis
## Get first event after 5 years of primary diagnosis
table(data$event == 1)
CA <- data[data$sndx == 4,]
CO <- data[data$sndx == 0,]
CA <- setDT(CA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

data <- rbind.data.frame(CA, CO)


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

alldata$start <- NA
alldata$end <- NA
### For those without event, start is agedx and end is Fu date === for SN, analysis starts from 5 years post DX (SNs within 5 years have been removed from the analysis)
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0] + 5
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx +5 and end is first event time; Achal: also added +5 in the line below
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


# SNs_py2=SNs_py
SNs_py2=SNs_py[SNs_py$evt1==1,]


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



## Add cubic spline of age attained
breaks = seq(5, 95, 22.5)
cp = quantile(SNs_py2$end, breaks/100, na.rm = T)
cs = cubic_spline(SNs_py2$end, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
SNs_py2 <- cbind.data.frame(SNs_py2, cs)

###############
## Model fit ##
###############
## 3. Sarcoma
fit_all <- glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                 sex +
                 cat_aa_class_dose +
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
results.Sarcoma <- data.frame(
  Coefficient = rownames(coefficients),
  Relative_Risk = round(relative_risks,1),
  CI_Lower = round(conf_int_low,1),
  CI_Upper = round(conf_int_high,1),
  P_Value = p_values
)

# Print the results
print(results.Sarcoma)


###################
## Breast cancer ##
###################
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

data <- data.all

# data$agelstcontact.original <- data$agelstcontact
data$gradedt <- data$evaldt ## SN grade date

## Age at SN diagnosis
data$AGE.ANY_SN <- time_length(interval(as.Date(data$dob), as.Date(data$gradedt)), "years")

## Attained age: age at last contact for cases is SN diagnosis date
data$agelstcontact[!is.na(data$evaldt)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

## define event (Any SN)
# 2= Meningioma, 4= Sarcoma, 5 = Breast cancer, 6 = Thyroid, 7 = NMSC
# data$event <- ifelse(data$sndx != 0, 1, 0) ## Any SN
data$event <- ifelse(data$sndx == 5, 1, 0) ## Breast cancer

## This is when we limit to first event
## Check how many within 5 years of primary cancer
data$AGE.ANY_SN.after.childhood.cancer.from.agedx <- data$AGE.ANY_SN - data$agedx
sum(data$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5, na.rm = T)
# 0 # Since it's 0, they all are after 5 years of primary diagnosis
## Get first event after 5 years of primary diagnosis
table(data$event == 1)
CA <- data[data$sndx == 5,]
CO <- data[data$sndx == 0,]
CA <- setDT(CA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

data <- rbind.data.frame(CA, CO)
dim(data)

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

alldata$start <- NA
alldata$end <- NA
### For those without event, start is agedx and end is Fu date === for SN, analysis starts from 5 years post DX (SNs within 5 years have been removed from the analysis)
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0] + 5
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx +5 and end is first event time; Achal: also added +5 in the line below
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
# SNs_py2$sex <- relevel(SNs_py2$sex, ref = "Male")
SNs_py2$braincat <- factor(SNs_py2$braincat)
SNs_py2$abdcat <- factor(SNs_py2$abdcat)
SNs_py2$chestcat <- factor(SNs_py2$chestcat)
SNs_py2$pelviscat <- factor(SNs_py2$pelviscat)
SNs_py2$neckcat <- factor(SNs_py2$neckcat)
SNs_py2$cat_anthra_jco_dose <- factor(SNs_py2$cat_anthra_jco_dose)
SNs_py2$cat_epitxn_dose <- factor(SNs_py2$cat_epitxn_dose)
SNs_py2$cat_aa_class_dose <- factor(SNs_py2$cat_aa_class_dose)



## Add cubic spline of age attained
breaks = seq(5, 95, 22.5)
cp = quantile(SNs_py2$end, breaks/100, na.rm = T)
cs = cubic_spline(SNs_py2$end, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
SNs_py2 <- cbind.data.frame(SNs_py2, cs)


# # library(splines)
# # # Define evenly spaced knots
# num_knots <- 5
# knots <- seq(min(SNs_py2$end), max(SNs_py2$end), length.out = num_knots + 2)[2:(num_knots + 1)]
# 
# # Create a cubic spline with 5 evenly spaced knots
# cubic_spline <- ns(SNs_py2$end, knots = knots)
# 
# # colnames(cubic_spline) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4", "AGE_AT_LAST_CONTACT.cs5")
# colnames(cubic_spline) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4", "AGE_AT_LAST_CONTACT.cs5", "AGE_AT_LAST_CONTACT.cs6")
# cs <- cubic_spline
# SNs_py2 <- cbind.data.frame(SNs_py2, cs)

###############
## Model fit ##
###############

# subset fameles only
SNs_py2 <- SNs_py2[SNs_py2$sex == "Female",]

## 4. Breast cancer
fit_all <- glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                 agedxcat +
                 chestcat.AN + cat_anthra_jco_dose.AN + 
                 EAS + AFR,
               family = "poisson", offset = log(SNs_py2$PY), data = SNs_py2)


# SNs_py2 <- SNs_py2 %>%
#   mutate(id = as.numeric(factor(sjlid)))
# # 
# SNs_py2 <- SNs_py2[c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4",
# "agedxcat", "chestcat", "cat_anthra_jco_dose", "event", "sjlid", "PY", "id")]
# # 
# SNs_py2 <- SNs_py2[complete.cases(SNs_py2), ]
# 
# fit_all <- geeglm(event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                     agedxcat + chestcat + cat_anthra_jco_dose, 
#                   family = "poisson", offset = log(SNs_py2$PY), id = id, corstr = "independence",  std.err = "san.se", data = SNs_py2)

summary_fit_all <- summary(fit_all)


# Extract coefficients, standard errors, and p-values from the model summary
coefficients <- coef(summary_fit_all)
std_errors <- coefficients[, "Std. Error"]
p_values <- coefficients[, "Pr(>|z|)"]

# Calculate relative risks (exponentiated coefficients)
relative_risks <- exp(coefficients[, "Estimate"])

# cbind.data.frame(rownames(coefficients), exp(coefficients[, "Estimate"]))

# Calculate 95% confidence intervals for relative risks
conf_int_low <- exp(coefficients[, "Estimate"] - 1.96 * std_errors)
conf_int_high <- exp(coefficients[, "Estimate"] + 1.96 * std_errors)

# Create a data frame to display the results
results.Breast <- data.frame(
  Coefficient = rownames(coefficients),
  Relative_Risk = round(relative_risks,1),
  CI_Lower = round(conf_int_low,1),
  CI_Upper = round(conf_int_high,1),
  P_Value = p_values
)

# Print the results
print(results.Breast)


#############
## Thyroid ##
#############
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

data <- data.all

# data$agelstcontact.original <- data$agelstcontact
data$gradedt <- data$evaldt ## SN grade date

## Age at SN diagnosis
data$AGE.ANY_SN <- time_length(interval(as.Date(data$dob), as.Date(data$gradedt)), "years")

## Attained age: age at last contact for cases is SN diagnosis date
data$agelstcontact[!is.na(data$evaldt)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

## define event (Any SN)
# 2= Meningioma, 4= Sarcoma, 5, Breast cancer, 6 = Thyroid, 7 = NMSC
data$event <- ifelse(data$sndx == 6, 1, 0) ## Sarcoma

## This is when we limit to first event
## Check how many within 5 years of primary cancer
data$AGE.ANY_SN.after.childhood.cancer.from.agedx <- data$AGE.ANY_SN - data$agedx
sum(data$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5, na.rm = T)
# 0 # Since it's 0, they all are after 5 years of primary diagnosis
## Get first event after 5 years of primary diagnosis
table(data$event == 1)
CA <- data[data$sndx == 6,]
CO <- data[data$sndx == 0,]
CA <- setDT(CA)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]

data <- rbind.data.frame(CA, CO)


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

alldata$start <- NA
alldata$end <- NA
### For those without event, start is agedx and end is Fu date === for SN, analysis starts from 5 years post DX (SNs within 5 years have been removed from the analysis)
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0] + 5
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx +5 and end is first event time; Achal: also added +5 in the line below
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


# SNs_py2=SNs_py
SNs_py2=SNs_py[SNs_py$evt1==1,]


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



## Add cubic spline of age attained
breaks = seq(5, 95, 22.5)
cp = quantile(SNs_py2$end, breaks/100, na.rm = T)
cs = cubic_spline(SNs_py2$end, knots = cp)

colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
SNs_py2 <- cbind.data.frame(SNs_py2, cs)

###############
## Model fit ##
###############
## 5. Thyroid
fit_all <- glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                 agedxcat + sex +
                 neckcat + cat_epitxn_dose,
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
results.Thyroid <- data.frame(
  Coefficient = rownames(coefficients),
  Relative_Risk = round(relative_risks,1),
  CI_Lower = round(conf_int_low,1),
  CI_Upper = round(conf_int_high,1),
  P_Value = p_values
)

# Print the results
print(results.Thyroid)

