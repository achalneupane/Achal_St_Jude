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




## MaxsegRT
PHENO.ANY_SN$maxsegrtdose.2020 <- radiation.2020$maxsegrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxsegrtdose.2018 <- radiation.2018$maxsegrtdose[match(PHENO.ANY_SN$sjlid, radiation.2018$sjlid)]
table(PHENO.ANY_SN$maxsegrtdose.2020 == PHENO.ANY_SN$maxsegrtdose) ## 2020
table(PHENO.ANY_SN$maxsegrtdose.2018 == PHENO.ANY_SN$maxsegrtdose) 


## Maxchest
PHENO.ANY_SN$maxchestrtdose.2020 <- radiation.2020$maxchestrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxchestrtdose.2018 <- radiation.2018$maxchestrtdose[match(PHENO.ANY_SN$sjlid, radiation.2018$sjlid)]
table(PHENO.ANY_SN$maxchestrtdose.2020 == PHENO.ANY_SN$maxchestrtdose) ## 2020
table(PHENO.ANY_SN$maxchestrtdose.2018 == PHENO.ANY_SN$maxchestrtdose) 


## Maxneck
PHENO.ANY_SN$maxneckrtdose.2020 <- radiation.2020$maxneckrtdose[match(PHENO.ANY_SN$sjlid, radiation.2020$sjlid)]
PHENO.ANY_SN$maxneckrtdose.2018 <- radiation.2018$maxneckrtdose[match(PHENO.ANY_SN$sjlid, radiation.2018$sjlid)]
table(PHENO.ANY_SN$maxneckrtdose.2020 == PHENO.ANY_SN$maxneckrtdose) ## 2020
table(PHENO.ANY_SN$maxneckrtdose.2018 == PHENO.ANY_SN$maxneckrtdose) 

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





##########################################################################################################
## 1.-----------------------------------maxsegRT
## Check rt and chemo in Qi's data
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]
dim(data.rtchemocheck)

data.rtchemocheck$braincat.AN <- cut(data.rtchemocheck$maxsegrtdose, breaks = c(0, 200, 1799, 2999, max(data.rtchemocheck$maxsegrtdose, na.rm = T)),
    labels = c(0:3),
    include.lowest = TRUE)
levels(data.rtchemocheck$braincat.AN) <- c(levels(data.rtchemocheck$braincat.AN), 4)
data.rtchemocheck$braincat.AN[is.na(data.rtchemocheck$braincat.AN)] <- 4

data.rtchemocheck$braincat.AN.2020 <- cut(data.rtchemocheck$maxsegrtdose.2020, breaks = c(0, 200, 1799, 2999, max(data.rtchemocheck$maxsegrtdose.2020, na.rm = T)),
                            labels = c(0:3),
                            include.lowest = TRUE)
levels(data.rtchemocheck$braincat.AN.2020) <- c(levels(data.rtchemocheck$braincat.AN.2020), 4)
data.rtchemocheck$braincat.AN.2020[is.na(data.rtchemocheck$braincat.AN.2020)] <- 4


data.rtchemocheck$braincat.AN.2018 <- cut(data.rtchemocheck$maxsegrtdose.2018, breaks = c(0, 200, 1799, 2999, max(data.rtchemocheck$maxsegrtdose.2018, na.rm = T)),
                                     labels = c(0:3),
                                     include.lowest = TRUE)
levels(data.rtchemocheck$braincat.AN.2018) <- c(levels(data.rtchemocheck$braincat.AN.2018), 4)
data.rtchemocheck$braincat.AN.2018[is.na(data.rtchemocheck$braincat.AN.2018)] <- 4

table(data.rtchemocheck$braincat)
# 0    1    2    3    4 
# 3218   18  556  279  331 

table(data.rtchemocheck$braincat.AN)
# 0    1    2    3    4 
# 2992   13  556  279  562

table(data.rtchemocheck$braincat.AN.2018)
# 0    1    2    3    4 
# 3258   18  568  509   49 

table(data.rtchemocheck$braincat.AN.2020)
# 0    1    2    3    4 
# 3255   18  568  509   52 

table(PHENO.ANY_SN$maxsegrtdose.category)
# None   >0-<18 >=18-<30     >=30  Unknown 
# 3255       18      568      508       52 

# table(data.rtchemocheck$maxsegrtdose.2018 == data.rtchemocheck$maxsegrtdose)
cc <- data.rtchemocheck[c("sjlid", "braincat", "braincat.AN", "braincat.AN.2020", "maxsegrtdose", "maxsegrtdose.2018", "maxsegrtdose.2020", "braindose")]

cc <- data[c("sjlid", "maxsegrtdose", "braindose", "braincat")]

## There are NAs but not labelled 4 in Qi's data
# 0    1    2    3    4 
# 3222   19  978  368  587
# Also, for example, IDs SJL4186407 and SJL5056617 are incorrectly categorized; and when there are missing values in SJL5244907 or SJL5245607..


##########################################################################################################
## 2.-----------------------------------maxneckRT
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]

data.rtchemocheck$neckcat.AN <- cut(data.rtchemocheck$neckmaxrtdose, breaks = c(0, 200, 1099, 1999, 2999, max(data.rtchemocheck$neckmaxrtdose, na.rm = T)),
    labels = c(0:4),
    include.lowest = TRUE)
levels(data.rtchemocheck$neckcat.AN) <- c(levels(data.rtchemocheck$neckcat.AN), 5)
data.rtchemocheck$neckcat.AN[is.na(data.rtchemocheck$neckcat.AN)] <- 5


data.rtchemocheck$neckcat.AN.2018 <- cut(data.rtchemocheck$maxneckrtdose.2018, breaks = c(0, 200, 1099, 1999, 2999, max(data.rtchemocheck$maxneckrtdose.2018, na.rm = T)),
                                         labels = c(0:4),
                                         include.lowest = TRUE)
levels(data.rtchemocheck$neckcat.AN.2018) <- c(levels(data.rtchemocheck$neckcat.AN.2018), 5)
data.rtchemocheck$neckcat.AN.2018[is.na(data.rtchemocheck$neckcat.AN.2018)] <- 5


data.rtchemocheck$neckcat.AN.2020 <- cut(data.rtchemocheck$maxneckrtdose.2020, breaks = c(0, 200, 1099, 1999, 2999, max(data.rtchemocheck$maxneckrtdose.2020, na.rm = T)),
                           labels = c(0:4),
                           include.lowest = TRUE)
levels(data.rtchemocheck$neckcat.AN.2020) <- c(levels(data.rtchemocheck$neckcat.AN.2020), 5)
data.rtchemocheck$neckcat.AN.2020[is.na(data.rtchemocheck$neckcat.AN.2020)] <- 5


table(data.rtchemocheck$neckcat)
# 0    1    2    3 
# 3397  233  368  404 

table(data.rtchemocheck$neckcat.AN)
# 0    1    2    3    4    5 
# 3247    5   64  306  223  557

table(data.rtchemocheck$neckcat.AN.2018)
# 0    1    2    3    4    5 
# 3565    6   90  426  270   45 

table(data.rtchemocheck$neckcat.AN.2020)
# 0    1    2    3    4    5 
# 3563    6   87  421  277   48 

table(PHENO.ANY_SN$maxneckrtdose.category)
# None   >0-<11 >=11-<20 >=20-<30     >=30  Unknown 
# 3562        6       87      421      277       48 

cc <- data.all[c("sjlid", "neckcat", "neckmaxrtdose", "maxneckrtdose.2018", "maxneckrtdose.2020")]

## Our data
table(PHENO.ANY_SN$maxneckrtdose.category)
# None   >0-<11 >=11-<20 >=20-<30     >=30  Unknown 
# 3561        6       90      426      270       48 


##########################################################################################################
## 3.-----------------------------------abdomenRT
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]


data.rtchemocheck$abdcat.AN <- cut(data.rtchemocheck$abdmaxrtdose, breaks = c(0, 200, 2999, max(data.rtchemocheck$abdmaxrtdose, na.rm = T)),
    labels = c(0:2),
    include.lowest = TRUE)
levels(data.rtchemocheck$abdcat.AN) <- c(levels(data.rtchemocheck$abdcat.AN), 3)
data.rtchemocheck$abdcat.AN [is.na(data.rtchemocheck$abdcat.AN)] <- 3


data.rtchemocheck$abdcat.AN.2018 <- cut(data.rtchemocheck$maxabdrtdose.2018, breaks = c(0, 200, 2999, max(data.rtchemocheck$maxabdrtdose.2018, na.rm = T)),
                          labels = c(0:2),
                          include.lowest = TRUE)
levels(data.rtchemocheck$abdcat.AN.2018) <- c(levels(data.rtchemocheck$abdcat.AN.2018), 3)
data.rtchemocheck$abdcat.AN.2018 [is.na(data.rtchemocheck$abdcat.AN.2018)] <- 3

data.rtchemocheck$abdcat.AN.2020 <- cut(data.rtchemocheck$maxabdrtdose.2020, breaks = c(0, 200, 2999, max(data.rtchemocheck$maxabdrtdose.2020, na.rm = T)),
                               labels = c(0:2),
                               include.lowest = TRUE)
levels(data.rtchemocheck$abdcat.AN.2020) <- c(levels(data.rtchemocheck$abdcat.AN.2020), 3)
data.rtchemocheck$abdcat.AN.2020 [is.na(data.rtchemocheck$abdcat.AN.2020)] <- 3

table(data.rtchemocheck$abdcat)
# 0    1    2    3 
# 3592  366  200  244 

table(data.rtchemocheck$abdcat.AN)
# 0    1    2    3 
# 3280  363  201  558 
table(data.rtchemocheck$abdcat.AN.2018)
# 0    1    2    3 
# 3572  522  263   45 
table(data.rtchemocheck$abdcat.AN.2020)
# 0    1    2    3 
# 3567  524  263   48 
## Our data
table(PHENO.ANY_SN$maxabdrtdose.category)
# None  >0-<30    >=30 Unknown 
# 3566     524     263      48 

##########################################################################################################
## 4.-----------------------------------chestRT
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]


data.rtchemocheck$chestcat.AN <- cut(data.rtchemocheck$chestmaxrtdose, breaks = c(0, 200, 1999, max(data.rtchemocheck$chestmaxrtdose, na.rm = T)),
                                   labels = c(0:2),
                                   include.lowest = TRUE)
levels(data.rtchemocheck$chestcat.AN) <- c(levels(data.rtchemocheck$chestcat.AN), 3)
data.rtchemocheck$chestcat.AN [is.na(data.rtchemocheck$chestcat.AN)] <- 3


data.rtchemocheck$chestcat.AN.2018 <- cut(data.rtchemocheck$maxchestrtdose.2018, breaks = c(0, 200, 1999, max(data.rtchemocheck$maxchestrtdose.2018, na.rm = T)),
                                        labels = c(0:2),
                                        include.lowest = TRUE)
levels(data.rtchemocheck$chestcat.AN.2018) <- c(levels(data.rtchemocheck$chestcat.AN.2018), 3)
data.rtchemocheck$chestcat.AN.2018 [is.na(data.rtchemocheck$chestcat.AN.2018)] <- 3

data.rtchemocheck$chestcat.AN.2020 <- cut(data.rtchemocheck$maxchestrtdose.2020, breaks = c(0, 200, 1999, max(data.rtchemocheck$maxchestrtdose.2020, na.rm = T)),
                                        labels = c(0:2),
                                        include.lowest = TRUE)
levels(data.rtchemocheck$chestcat.AN.2020) <- c(levels(data.rtchemocheck$chestcat.AN.2020), 3)
data.rtchemocheck$chestcat.AN.2020 [is.na(data.rtchemocheck$chestcat.AN.2020)] <- 3

table(data.rtchemocheck$chestcat)
# 0    1    2    3 
# 3463  109  518  312 

table(data.rtchemocheck$chestcat.AN)
# 0    1    2    3 
# 3221  105  518  558 
table(data.rtchemocheck$chestcat.AN.2018)
# 0    1    2    3 
# 3502  149  705   46 
table(data.rtchemocheck$chestcat.AN.2020)
# 0    1    2    3 
# 3501  149  703   49 
## Our data
table(PHENO.ANY_SN$maxchestrtdose.category)
# None  >0-<20    >=20 Unknown 
# 3500     149     703      49 

table(PHENO.ANY_SN$maxchestrtdose == PHENO.ANY_SN$maxchestrtdose.2020)
# table(data.rtchemocheck$chestmaxrtdose == data.rtchemocheck$maxchestrtdose.2020)
table(!is.na(data.rtchemocheck$chestmaxrtdose) & data.rtchemocheck$chestmaxrtdose == data.rtchemocheck$maxchestrtdose.2020)

cc <- data.rtchemocheck[c("sjlid", "chestmaxrtdose", "chestdose", "maxchestrtdose.2018", "maxchestrtdose.2020", "chestcat", "chestcat.AN", "chestcat.AN.2018", "chestcat.AN.2020")]
## Take screenshot of cc

cc <- data[c("sjlid", "chestmaxrtdose", "chestdose", "chestcat")]

##########################################################################################################
## 5.-----------------------------------pelvisRT
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]


data.rtchemocheck$pelviscat.AN <- cut(data.rtchemocheck$pelvismaxrtdose, breaks = c(0, 200, 1999, max(data.rtchemocheck$pelvismaxrtdose, na.rm = T)),
                                     labels = c(0:2),
                                     include.lowest = TRUE)
levels(data.rtchemocheck$pelviscat.AN) <- c(levels(data.rtchemocheck$pelviscat.AN), 3)
data.rtchemocheck$pelviscat.AN [is.na(data.rtchemocheck$pelviscat.AN)] <- 3


data.rtchemocheck$pelviscat.AN.2018 <- cut(data.rtchemocheck$maxpelvisrtdose.2018, breaks = c(0, 200, 1999, max(data.rtchemocheck$maxpelvisrtdose.2018, na.rm = T)),
                                          labels = c(0:2),
                                          include.lowest = TRUE)
levels(data.rtchemocheck$pelviscat.AN.2018) <- c(levels(data.rtchemocheck$pelviscat.AN.2018), 3)
data.rtchemocheck$pelviscat.AN.2018 [is.na(data.rtchemocheck$pelviscat.AN.2018)] <- 3

data.rtchemocheck$pelviscat.AN.2020 <- cut(data.rtchemocheck$maxpelvisrtdose.2020, breaks = c(0, 200, 1999, max(data.rtchemocheck$maxpelvisrtdose.2020, na.rm = T)),
                                          labels = c(0:2),
                                          include.lowest = TRUE)
levels(data.rtchemocheck$pelviscat.AN.2020) <- c(levels(data.rtchemocheck$pelviscat.AN.2020), 3)
data.rtchemocheck$pelviscat.AN.2020 [is.na(data.rtchemocheck$pelviscat.AN.2020)] <- 3

table(data.rtchemocheck$pelviscat)
# 0    1    2    3 
# 3685  129  377  211 

table(data.rtchemocheck$pelviscat.AN)
# 0    1    2    3 
# 3341  126  377  558 
table(data.rtchemocheck$pelviscat.AN.2018)
# 0    1    2    3 
# 3673  169  518   42 
table(data.rtchemocheck$pelviscat.AN.2020)
# 0    1    2    3 
# 3669  169  519   45 
## Our data
table(PHENO.ANY_SN$maxpelvisrtdose.category)
# None  >0-<20    >=20 Unknown 
# 3668     169     519      45 

table(PHENO.ANY_SN$maxpelvisrtdose == PHENO.ANY_SN$maxpelvisrtdose.2020)
# table(data.rtchemocheck$pelvismaxrtdose == data.rtchemocheck$maxpelvisrtdose.2020)
table(!is.na(data.rtchemocheck$pelvismaxrtdose) & data.rtchemocheck$pelvismaxrtdose == data.rtchemocheck$maxpelvisrtdose.2020)

cc <- data.rtchemocheck[c("sjlid", "pelvismaxrtdose", "maxpelvisrtdose.2018", "maxpelvisrtdose.2020", "pelviscat", "pelviscat.AN", "pelviscat.AN.2018", "pelviscat.AN.2020")]
## Take screenshot of cc

##########################################################################################################
## 6.-----------------------------------aa_class_dose
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]


TERT = unname(quantile(data.rtchemocheck$aa_class_dose[data.rtchemocheck$aa_class_dose !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_aa_class_dose.AN <- cut(data.rtchemocheck$aa_class_dose, breaks = c(-Inf, 0, TERT),
                                               labels = c("None", "1st", "2nd", "3rd"),
                                               include.lowest = TRUE)



TERT = unname(quantile(data.rtchemocheck$aa_class_dose.2018[data.rtchemocheck$aa_class_dose.2018 !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_aa_class_dose.AN.2018 <- cut(data.rtchemocheck$aa_class_dose.2018, breaks = c(-Inf, 0, TERT),
                                                   labels = c("None", "1st", "2nd", "3rd"),
                                                   include.lowest = TRUE)

TERT = unname(quantile(data.rtchemocheck$aa_class_dose.2020[data.rtchemocheck$aa_class_dose.2020 !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_aa_class_dose.AN.2020 <- cut(data.rtchemocheck$aa_class_dose.2020, breaks = c(-Inf, 0, TERT),
                                                   labels = c("None", "1st", "2nd", "3rd"),
                                                   include.lowest = TRUE)


table(data.rtchemocheck$cat_aa_class_dose)
# 0    1    2    3 
# 1911  820  824  820
table(data.rtchemocheck$cat_aa_class_dose.AN)
# None  1st  2nd  3rd 
# 1911  822  821  821 
table(data.rtchemocheck$cat_aa_class_dose.AN.2018)
# None  1st  2nd  3rd 
# 1920  823  823  823
table(data.rtchemocheck$cat_aa_class_dose.AN.2020)
# None  1st  2nd  3rd 
# 1943  819  812  816 
table(PHENO.ANY_SN$aa_class_dose_5.category)
# None  1st  2nd  3rd 
# 1943  819  812  815


data.rtchemocheck$aa_class_dose.AN.2020 <- PHENO.ANY_SN$aa_class_dose_5.2020[match(data.rtchemocheck$sjlid, PHENO.ANY_SN$sjlid)]
data.rtchemocheck$aa_class_dose.AN.2018 <- PHENO.ANY_SN$aa_class_dose_5.2018[match(data.rtchemocheck$sjlid, PHENO.ANY_SN$sjlid)]
table(data.rtchemocheck$aa_class_dose.AN.2018 == data.rtchemocheck$aa_class_dose)
table(data.rtchemocheck$aa_class_dose.AN.2020 == data.rtchemocheck$aa_class_dose)

table(data.rtchemocheck$aa_class_dose_5.category)
# None  1st  2nd  3rd 
# 1943  819  812  815

cc <- data.rtchemocheck[c("sjlid", "aa_class_dose", "aa_class_dose.AN.2018", "aa_class_dose.AN.2020", "cat_aa_class_dose", "cat_aa_class_dose.AN", "cat_aa_class_dose.AN.2018", "cat_aa_class_dose.AN.2020")]


##########################################################################################################
## 7.-----------------------------------anthra_jco_dose
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]


TERT = unname(quantile(data.rtchemocheck$anthra_jco_dose[data.rtchemocheck$anthra_jco_dose !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_anthra_jco_dose.AN <- cut(data.rtchemocheck$anthra_jco_dose, breaks = c(-Inf, 0, TERT),
                                              labels = c("None", "1st", "2nd", "3rd"),
                                              include.lowest = TRUE)


TERT = unname(quantile(data.rtchemocheck$anthra_jco_dose.2018[data.rtchemocheck$anthra_jco_dose.2018 !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_anthra_jco_dose.AN.2018 <- cut(data.rtchemocheck$anthra_jco_dose.2018, breaks = c(-Inf, 0, TERT),
                                                   labels = c("None", "1st", "2nd", "3rd"),
                                                   include.lowest = TRUE)

TERT = unname(quantile(data.rtchemocheck$anthra_jco_dose.2020[data.rtchemocheck$anthra_jco_dose.2020 !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_anthra_jco_dose.AN.2020 <- cut(data.rtchemocheck$anthra_jco_dose.2020, breaks = c(-Inf, 0, TERT),
                                                   labels = c("None", "1st", "2nd", "3rd"),
                                                   include.lowest = TRUE)


table(data.rtchemocheck$cat_anthra_jco_dose)
# 0    1    2    3 
# 1938  815  815  815
table(data.rtchemocheck$cat_anthra_jco_dose.AN)
# None  1st  2nd  3rd 
# 1938  815  815  815

table(data.rtchemocheck$cat_anthra_jco_dose.AN.2018)
# None  1st  2nd  3rd 
# 1947  817  816  817
table(data.rtchemocheck$cat_anthra_jco_dose.AN.2020)
# None  1st  2nd  3rd 
# 1973  809  808  808
table(PHENO.ANY_SN$anthra_jco_dose_5.category)
# None  1st  2nd  3rd 
# 1972  809  808  808 


data.rtchemocheck$anthra_jco_dose.AN.2020 <- PHENO.ANY_SN$anthra_jco_dose_5.2020[match(data.rtchemocheck$sjlid, PHENO.ANY_SN$sjlid)]
data.rtchemocheck$anthra_jco_dose.AN.2018 <- PHENO.ANY_SN$anthra_jco_dose_5.2018[match(data.rtchemocheck$sjlid, PHENO.ANY_SN$sjlid)]
table(data.rtchemocheck$anthra_jco_dose.AN.2018 == data.rtchemocheck$anthra_jco_dose)
table(data.rtchemocheck$anthra_jco_dose.AN.2020 == data.rtchemocheck$anthra_jco_dose)

cc <- data.rtchemocheck[c("sjlid", "anthra_jco_dose", "anthra_jco_dose.AN.2018", "anthra_jco_dose.AN.2020", "cat_anthra_jco_dose", "cat_anthra_jco_dose.AN", "cat_anthra_jco_dose.AN.2018", "cat_anthra_jco_dose.AN.2020")]


##########################################################################################################
## 8.-----------------------------------epitxn_dose
data.rtchemocheck <- data.all[!duplicated(data.all$sjlid),]


TERT = unname(quantile(data.rtchemocheck$epitxn_dose[data.rtchemocheck$epitxn_dose !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_epitxn_dose.AN <- cut(data.rtchemocheck$epitxn_dose, breaks = c(-Inf, 0, TERT),
                                                labels = c("None", "1st", "2nd", "3rd"),
                                                include.lowest = TRUE)


TERT = unname(quantile(data.rtchemocheck$epitxn_dose.2018[data.rtchemocheck$epitxn_dose.2018 !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_epitxn_dose.AN.2018 <- cut(data.rtchemocheck$epitxn_dose.2018, breaks = c(-Inf, 0, TERT),
                                                     labels = c("None", "1st", "2nd", "3rd"),
                                                     include.lowest = TRUE)

TERT = unname(quantile(data.rtchemocheck$epitxn_dose.2020[data.rtchemocheck$epitxn_dose.2020 !=0], c(1/3, 2/3, 1), na.rm = T))
data.rtchemocheck$cat_epitxn_dose.AN.2020 <- cut(data.rtchemocheck$epitxn_dose.2020, breaks = c(-Inf, 0, TERT),
                                                     labels = c("None", "1st", "2nd", "3rd"),
                                                     include.lowest = TRUE)


table(data.rtchemocheck$cat_epitxn_dose)
# 0    1    2    3 
# 2896  494  496  494 
table(data.rtchemocheck$cat_epitxn_dose.AN)
# 0    1    2    3 
# 2896  495  494  495
table(data.rtchemocheck$cat_epitxn_dose.AN.2018)
# None  1st  2nd  3rd 
# 2906  496  496  496 
table(data.rtchemocheck$cat_epitxn_dose.AN.2020)
# None  1st  2nd  3rd 
# 2906  496  495  496 
table(PHENO.ANY_SN$epitxn_dose_5.category)
# None  1st  2nd  3rd 
# 2906  496  495  495 


data.rtchemocheck$epitxn_dose.AN.2020 <- PHENO.ANY_SN$epitxn_dose_5.2020[match(data.rtchemocheck$sjlid, PHENO.ANY_SN$sjlid)]
data.rtchemocheck$epitxn_dose.AN.2018 <- PHENO.ANY_SN$epitxn_dose_5.2018[match(data.rtchemocheck$sjlid, PHENO.ANY_SN$sjlid)]
table(data.rtchemocheck$epitxn_dose.AN.2018 == data.rtchemocheck$epitxn_dose)
table(data.rtchemocheck$epitxn_dose.AN.2020 == data.rtchemocheck$epitxn_dose)

cc <- data.rtchemocheck[c("sjlid", "epitxn_dose", "epitxn_dose.AN.2018", "epitxn_dose.AN.2020", "cat_epitxn_dose", "cat_epitxn_dose.AN", "cat_epitxn_dose.AN.2018", "cat_epitxn_dose.AN.2020")]




data.rtchemocheck$agedx
data.rtchemocheck$agedxcat.AN <- data.rtchemocheck$agedx
data.rtchemocheck$agedxcat.AN[data.rtchemocheck$agedx >= 0 & data.rtchemocheck$agedx < 5 ] <- "0-4"
data.rtchemocheck$agedxcat.AN[data.rtchemocheck$agedx >= 5 & data.rtchemocheck$agedx < 10 ] <- "5-9"
data.rtchemocheck$agedxcat.AN[data.rtchemocheck$agedx >= 10 & data.rtchemocheck$agedx < 15 ] <- "10-14"
data.rtchemocheck$agedxcat.AN[data.rtchemocheck$agedx >= 15 ] <- ">=15"
data.rtchemocheck$agedxcat.AN <- factor(data.rtchemocheck$agedxcat.AN, levels = c("0-4", "5-9", "10-14", ">=15")) # first level will be treated as reference
table(data.rtchemocheck$agedxcat.AN)
table(data.rtchemocheck$agedxcat)
table(PHENO.ANY_SN$AGE_AT_DIAGNOSIS)

table(data.rtchemocheck$sex)
table(PHENO.ANY_SN$gender)
