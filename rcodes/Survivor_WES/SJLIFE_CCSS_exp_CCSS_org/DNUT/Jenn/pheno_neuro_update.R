library(haven)
library(dplyr)
library(tidyverse)
library(readxl)

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Neurocognition/")

# read in neurocognitive outcomes
neuro <- read.csv("neuro_all_data.csv", header = TRUE,sep = ',')
neuro$doe <- as.Date(neuro$doe, "%m/%d/%Y")

# sort by age at exam - pick baseline for each individual
ntest3 <- neuro %>% group_by(sjlid) 

ncounts <- ntest3 %>% summarise(n=n())

ncounts$n <- as.numeric(ncounts$n)

 # loop to read line by line to see if there is one or multiple measures for surgery
col_names <- c("sjlid","doe","ageexam","CPTOmZs","CPTVarZs","DigitFwdZs","DigitFwdSpanZs","TrailAZs","CVLTSDFRZs","CVLTLDFRZs",
               "CVLTTotZs","VisualSelZs","TrailBZs","VerFluZs","DigitBkwdZs","DigitBkwdSpanZs","CPTComZs","CodingZs","SymSrchZs",
               "GPDDHZs","WASIVocZs","WASIMatZs","WASIFSIQZs","WASIFSIQ","WJLWIDZs","WJCalcZs","ReyCopyZS","current_age")
base_neuro <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
colnames(base_neuro) <- col_names
base_neuro$sjlid <- as.character(base_neuro$sjlid)
base_neuro$doe <- as.Date(base_neuro$doe)
base_neuro[,3:28] <- lapply(base_neuro[,3:28], as.numeric)


recent_neuro <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
colnames(recent_neuro) <- col_names
recent_neuro$sjlid <- as.character(recent_neuro$sjlid)
recent_neuro$doe <- as.Date(recent_neuro$doe)
recent_neuro[,3:28] <- lapply(recent_neuro[,3:28], as.numeric)


for (i in 1:length(ncounts$n)) {
  # get line number to print subid
  cat(i)
  # id to grep in test3 file
  sid <- ncounts$sjlid[i]
    
  # ungroup test3 so I can add lines to output
  ntest3 <- ntest3 %>% ungroup()
  filtered_row <- grep(sid,ntest3$sjlid,value=FALSE)
  sub_lines <- neuro[filtered_row,]
  sub_lines <- sub_lines[order(sub_lines$ageexam),]
  
  # create file with matching line from test3
  base_neuro <- base_neuro %>% add_row(sub_lines[1,])
  
  
  # now sort in descending order to get most recent measure
  sub_lines <- sub_lines[order(sub_lines$ageexam, decreasing = TRUE),]
  
  # create file with matching line from test3
  recent_neuro <- recent_neuro %>% add_row(sub_lines[1,])
}

base_neuro$doe <- as.Date(base_neuro$doe)
recent_neuro$doe <- as.Date(recent_neuro$doe)





# chemotherapy variables
# path for mac: /Volumes/sjcommon/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat
chemo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat")
#chemo2 <- select(chemo,sjlid,methotrexate_HD_dose_5,methotrexate_HD_dose_10,methotrexate_HD_dose_any,methotrexate_HD_dose_prim,
#                 methotrexate_it_io_dose_5,methotrexate_it_io_dose_10,methotrexate_it_io_dose_any,methotrexate_it_io_dose_prim,
#                 cytarabine_HD_dose_5,cytarabine_HD_dose_10,cytarabine_HD_dose_any,cytarabine_HD_dose_prim, cytarabine_dose_5,
#                 cytarabine_dose_10,cytarabine_dose_any,cytarabine_dose_prim,alkylating_dose_5, alkylating_dose_10,
#                 alkylating_dose_prim,alkylating_dose_any,dexamethasone_dose_5,dexamethasone_dose_10,dexamethasone_dose_any,
#                 dexamethasone_dose_prim, anthracyclines_dose_5,anthracyclines_dose_10,anthracyclines_dose_any,
#                 anthracyclines_dose_prim)

chemo2 <- select(chemo,sjlid,methotrexate_HD_dose_5,methotrexate_it_io_dose_5,cytarabine_HD_dose_5,cytarabine_dose_5,
                 alkylating_dose_5,dexamethasone_dose_5,anthracyclines_dose_5)

######################################try later - taking a long time
# milli laboratory measures 
lab <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/milli_labs.sas7bdat")
lab2 <- select(lab,sjlid,labdt,cortisol)


# sort by age at exam - pick baseline for each individual
ntest3 <- lab2 %>% group_by(sjlid) 

ncounts <- ntest3 %>% summarise(n=n())

ncounts$n <- as.numeric(ncounts$n)

# select baseline
col_names <- c("sjlid","labdt","cortisol","diff")
base_lab <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
colnames(base_lab) <- col_names
base_lab$sjlid <- as.character(base_lab$sjlid)
base_lab$labdt <- as.Date(base_lab$labdt)
base_lab$cortisol <- as.character(base_lab$cortisol)
base_lab$diff <- as.character(base_lab$diff)
base_lab$diff <- as.difftime(base_lab$diff)


for (i in 1:length(ncounts$n)) {
  # get line number to print subid
  cat(i)
  # id to grep in test3 file
  sid <- ncounts$sjlid[i]
  
  # ungroup test3 so I can add lines to output
  ntest3 <- ntest3 %>% ungroup()
  filtered_row <- grep(sid,ntest3$sjlid,value=FALSE)
  neuro_base <- grep(sid,base_neuro$sjlid,value = FALSE)
  neuro <- base_neuro[neuro_base,]
  sub_lines <- lab2[filtered_row,]
  sub_lines$diff = abs(difftime(neuro$doe[1],sub_lines$labdt,units = "weeks"))
  
  # remove any lines with missing cortisol values
  sub_lines <- sub_lines[sub_lines$cortisol != "",]
  sub_lines <- sub_lines %>% arrange(diff)

  if (nrow(sub_lines) >= 1) {
  # create file with matching line from test3
  base_lab <- base_lab %>% add_row(sub_lines[1,])
  } else {
  base_lab <- base_lab %>% add_row(sjlid = sid, labdt = NA, cortisol = NA, diff = NA)
  }
}

base_lab <- select(base_lab,sjlid,cortisol)


# radiation measures ----- should maybe use head from radiation_dose?
rt <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiationsum_yn.sas7bdat")
rt2 <- select(rt,sjlid,Cranial_5)

# surgery data
surg <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/surgery.sas7bdat")
surg2 <- select(surg,sjlid,surgdt,brainsurg,brainsurgccss,ICD9_1,ICD9_2)


# add VP_shunt variable: 0 = no, 1 = yes
surg2$VPshunt <- ifelse(surg2$ICD9_1 == "02.34" | surg2$ICD9_2 == "02.34" | surg2$ICD9_1 == "02.39" | surg2$ICD9_2 == "02.39",1,0)

surg2$VPshunt[surg2$ICD9_1 == "02.34"] <- 1


# test loop for reading line by line
# make variable for id
test3 <- surg2 %>% group_by(sjlid) 

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


# loop to read line by line to see if there is one or multiple measures for surgery
for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  if (counts$n[i] == 1) {
    # get line number to print subid
    cat(i)
    # id to grep in test3 file
    sid <- counts$sjlid[i]
    
    # ungroup test3 so I can add lines to output
    test3 <- test3 %>% ungroup()
    filtered_row <- grep(sid,test3$sjlid,value=FALSE)
    # create file with matching line from test3
    done <- test3[filtered_row,]
    if (i == 1) { 
      outputsg <- done
    } else {
      done <- done %>% ungroup()
      outputsg <- outputsg %>% add_row(done)
    }
  } else {
    # has more than one measure of surgery
    cat(i)
    sid <- counts$sjlid[i]
    filtered <- grep(sid,test3$sjlid, value=FALSE)
    matches <- test3[filtered,]
    
    ones <- matches[matches$brainsurgccss == 1,]
    
    if( nrow(ones) > 0) {
      stat <- ones %>% arrange(surgdt)
      ones <- ones %>% ungroup()
      if ( i == 1) {
        outputsg <- ones[1,]
      } else {
        test3 <- test3 %>% ungroup()
        outputsg <- outputsg %>% add_row(ones[1,])
      }
    } else {
      zeros <- matches[matches$brainsurgccss == 0,]
      zeros <- zeros %>% arrange(surgdt)
      zeros <- zeros %>% ungroup()
      if ( i == 1) {
        outputsg <- zeros[1,]
      } else {
        outputsg <- outputsg %>% add_row(zeros[1,])
      }
    }
  }
}

outputsg <- select(outputsg,sjlid,brainsurgccss)

# Vp shunt
test3 <- surg2 %>% group_by(sjlid) 

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


# loop to read line by line to see if there is one or multiple measures for shunt
for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  if (counts$n[i] == 1) {
    # get line number to print subid
    cat(i)
    # id to grep in test3 file
    sid <- counts$sjlid[i]
    
    # ungroup test3 so I can add lines to output
    test3 <- test3 %>% ungroup()
    filtered_row <- grep(sid,test3$sjlid,value=FALSE)
    # create file with matching line from test3
    done <- test3[filtered_row,]
    if (i == 1) { 
      outputvp <- done
    } else {
      done <- done %>% ungroup()
      outputvp <- outputvp %>% add_row(done)
    }
  } else {
    # has more than one measure of surgery
    cat(i)
    sid <- counts$sjlid[i]
    filtered <- grep(sid,test3$sjlid, value=FALSE)
    matches <- test3[filtered,]
    
    ones <- matches[matches$VPshunt == 1,]
    
    if( nrow(ones) > 0) {
      stat <- ones %>% arrange(surgdt)
      ones <- ones %>% ungroup()
      if ( i == 1) {
        outputvp <- ones[1,]
      } else {
        test3 <- test3 %>% ungroup()
        outputvp <- outputvp %>% add_row(ones[1,])
      }
    } else {
      zeros <- matches[matches$VPshunt == 0,]
      zeros <- zeros %>% arrange(surgdt)
      zeros <- zeros %>% ungroup()
      if ( i == 1) {
        outputvp <- zeros[1,]
      } else {
        outputvp <- outputvp %>% add_row(zeros[1,])
      }
    }
  }
}

outputvp <- select(outputvp,sjlid,VPshunt)

surg_out <- merge(outputsg,outputvp)

surg_out$Surg_Shunt <- ifelse(surg_out$brainsurgccss == 1 | surg_out$brainsurgccss == 1,1,0)

surg_out <- select(surg_out,sjlid,Surg_Shunt)





# demographics
demo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
demo2 <- select(demo,sjlid,racegrp2,gender)

demo2$sex <- NA
demo2$sex[demo2$gender == "Female"] <- 2
demo2$sex[demo2$gender == "Male"] <- 1

demo2 <- select(demo2,sjlid,sex,racegrp2)
demo2$racegrp2 <- sub(",","",demo2$racegrp2)
demo2$racegrp2 <- gsub(" ","_",demo2$racegrp2)


# diagnosis
diag <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat")
diag2 <- select(diag,sjlid,agedx,primdx,diagsite,solidtumorsite)

diag3 <- diag2[diag2$primdx == 1 & !is.na(diag2$primdx),]
diag3$diagsite <- gsub(",","",diag3$diagsite)
diag3$solidtumorsite <- gsub(",","",diag3$solidtumorsite)

# CTCAE ------- come back to tomorrow - MACE events
grade <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")

cardio <- grade[grade$organsys == "Cardiovascular",]

# select worst grade
test3 <- cardio %>% group_by(sjlid) 

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  if (counts$n[i] == 1) {
    # get line number to print subid
    cat(i)
    # id to grep in test3 file
    sid <- counts$sjlid[i]
    
    # ungroup test3 so I can add lines to output
    test3 <- test3 %>% ungroup()
    filtered_row <- grep(sid,test3$sjlid,value=FALSE)
    # create file with matching line from test3
    done <- test3[filtered_row,]
    if (i == 1) { 
      out_cardio <- done
    } else {
      done <- done %>% ungroup()
      out_cardio <- out_cardio %>% add_row(done)
    }
  } else {
    # has more than one measure of surgery
    cat(i)
    sid <- counts$sjlid[i]
    filtered <- grep(sid,test3$sjlid, value=FALSE)
    matches <- test3[filtered,]
    matches <- matches %>% arrange(desc(matches$grade))
    matches <- matches %>% ungroup()
    if ( i == 1) {
      out_cardio <- matches[1,]
    } else {
      test3 <- test3 %>% ungroup()
      out_cardio <- out_cardio %>% add_row(matches[1,])
    }
  }
}

out_cardio$cardio_grade <- out_cardio$grade
out_cardio$cardio_CHC <- ifelse(out_cardio$grade >= 2, 1,0)

cardio_CHC <- select(out_cardio,sjlid,cardio_CHC)


# pulmonary
pulm <- grade[grade$organsys == "Pulmonary",]

# select worst grade

test3 <- pulm %>% group_by(sjlid) 

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  if (counts$n[i] == 1) {
    # get line number to print subid
    cat(i)
    # id to grep in test3 file
    sid <- counts$sjlid[i]
    
    # ungroup test3 so I can add lines to output
    test3 <- test3 %>% ungroup()
    filtered_row <- grep(sid,test3$sjlid,value=FALSE)
    # create file with matching line from test3
    done <- test3[filtered_row,]
    if (i == 1) { 
      out_pulm <- done
    } else {
      done <- done %>% ungroup()
      out_pulm <- out_pulm %>% add_row(done)
    }
  } else {
    # has more than one measure of surgery
    cat(i)
    sid <- counts$sjlid[i]
    filtered <- grep(sid,test3$sjlid, value=FALSE)
    matches <- test3[filtered,]
    matches <- matches %>% arrange(desc(matches$grade))
    matches <- matches %>% ungroup()
    if ( i == 1) {
      out_pulm <- matches[1,]
    } else {
      test3 <- test3 %>% ungroup()
      out_pulm <- out_pulm %>% add_row(matches[1,])
    }
  }
}

out_pulm$pulm_grade <- out_pulm$grade
out_pulm$pulm_CHC <- ifelse(out_pulm$grade >= 2, 1,0)

pulm_CHC <- select(out_pulm,sjlid,pulm_grade,pulm_CHC)


# endocrine
endo <- grade[grade$organsys == "Endocrine Diagnosis" | grade$organsys == "Endocrine Measure",]
# select worst grade

test3 <- endo %>% group_by(sjlid) 

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  if (counts$n[i] == 1) {
    # get line number to print subid
    cat(i)
    # id to grep in test3 file
    sid <- counts$sjlid[i]
    
    # ungroup test3 so I can add lines to output
    test3 <- test3 %>% ungroup()
    filtered_row <- grep(sid,test3$sjlid,value=FALSE)
    # create file with matching line from test3
    done <- test3[filtered_row,]
    if (i == 1) { 
      out_endo <- done
    } else {
      done <- done %>% ungroup()
      out_endo <- out_endo %>% add_row(done)
    }
  } else {
    # has more than one measure of surgery
    cat(i)
    sid <- counts$sjlid[i]
    filtered <- grep(sid,test3$sjlid, value=FALSE)
    matches <- test3[filtered,]
    matches <- matches %>% arrange(desc(matches$grade))
    matches <- matches %>% ungroup()
    if ( i == 1) {
      out_endo <- matches[1,]
    } else {
      test3 <- test3 %>% ungroup()
      out_endo <- out_endo %>% add_row(matches[1,])
    }
  }
}

out_endo$endo_grade <- out_endo$grade
out_endo$endo_CHC <- ifelse(out_endo$grade >= 2, 1,0)

endo_CHC <- select(out_endo,sjlid,endo_grade,endo_CHC)

# neuro - CHC
neur <- grade[grade$organsys == "Neurology",]

# select worst grade

test3 <- neur %>% group_by(sjlid) 

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  if (counts$n[i] == 1) {
    # get line number to print subid
    cat(i)
    # id to grep in test3 file
    sid <- counts$sjlid[i]
    
    # ungroup test3 so I can add lines to output
    test3 <- test3 %>% ungroup()
    filtered_row <- grep(sid,test3$sjlid,value=FALSE)
    # create file with matching line from test3
    done <- test3[filtered_row,]
    if (i == 1) { 
      out_neuro <- done
    } else {
      done <- done %>% ungroup()
      out_neuro <- out_neuro %>% add_row(done)
    }
  } else {
    # has more than one measure of surgery
    cat(i)
    sid <- counts$sjlid[i]
    filtered <- grep(sid,test3$sjlid, value=FALSE)
    matches <- test3[filtered,]
    matches <- matches %>% arrange(desc(matches$grade))
    matches <- matches %>% ungroup()
    if ( i == 1) {
      out_neuro <- matches[1,]
    } else {
      test3 <- test3 %>% ungroup()
      out_neuro <- out_neuro %>% add_row(matches[1,])
    }
  }
}

out_neuro$neuro_grade <- out_neuro$grade
out_neuro$neuro_CHC <- ifelse(out_neuro$grade >= 2, 1,0)

neuro_CHC <- select(out_neuro,sjlid,neuro_grade,neuro_CHC)

merge1 <- merge(cardio_CHC,pulm_CHC, all = TRUE)
merge2 <- merge(merge1,endo_CHC, all = TRUE)
ctcae <- merge(merge2,neuro_CHC, all = TRUE)

# chemo yn - skipped - check if needed after
chemoyn <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_yn.sas7bdat")
chemoyn2 <- select(chemoyn,sjlid,anthracyclines_5)


# smoking stat - need to recode as current, former, never
smk <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Survey Data/adult_healthhabits.sas7bdat")
smk2 <- select(smk,SJLIFEID,datecomp,smnow,evsm)

smk2$smk_stat <- NULL

smk2$smk_stat[smk2$evsm == 2] <- "Never"
smk2$smk_stat[smk2$evsm == 1 & smk2$smnow == 2] <- "Former"
smk2$smk_stat[smk2$evsm == 1 & smk2$smnow == 1] <- "Current"

smk <- smk2[!is.na(smk2$smk_stat),]
colnames(smk) <- c("sjlid","datecomp","smnow","evsm","smk_stat")

# select smoking data closest to baseline neuro
test3 <- smk %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$datecomp)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_base <- grep(sid,base_neuro$sjlid, value = FALSE)
  neuro <- base_neuro[neuro_base,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    smk_out <- matches[1,]
  } else {
    smk_out <- smk_out %>% add_row(matches[1,])
  }
}

smk_base <- select(smk_out,sjlid,smk_stat)

# smoking closest to recent neuro
for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_recent <- grep(sid,recent_neuro$sjlid, value = FALSE)
  neuro <- recent_neuro[neuro_recent,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  
  if ( i == 1) {
    smk_rec <- matches[1,]
  } else {
    smk_rec <- smk_rec %>% add_row(matches[1,])
  }
}

smk_recent <- select(smk_rec,sjlid,smk_stat)




# adolescent behavior

kbeh <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Survey Data/adolescent_behavior.sas7bdat")
kbeh2 <- select(kbeh,SJLIFEID,datecomp,grade,anxdx,depdx)

#kbeh3 <- kbeh2[!is.na(kbeh2$grade) | !is.na(kbeh2$anxdx) | !is.na(kbeh2$depdx),]


# adult behavior
abeh <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Survey Data/adult_behavior.sas7bdat")
abeh2 <- select(abeh,SJLIFEID,datecomp,grade,anxdx,depdx)

behav <- rbind(kbeh2,abeh2)
colnames(behav) <- c("sjlid","datecomp","grade","anxdx","depdx")

# select behavior measurements
grd <- behav[!is.na(behav$grade),]

# select grade data closest to baseline neuro
test3 <- grd %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$datecomp)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_base <- grep(sid,base_neuro$sjlid, value = FALSE)
  neuro <- base_neuro[neuro_base,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    grd_out <- matches[1,]
  } else {
    grd_out <- grd_out %>% add_row(matches[1,])
  }
}

grd_base <- select(grd_out,sjlid,grade)

# anxiety
anx <- behav[!is.na(behav$anxdx),]

# select anxiety data closest to baseline neuro
test3 <- anx %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$datecomp)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_base <- grep(sid,base_neuro$sjlid, value = FALSE)
  neuro <- base_neuro[neuro_base,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    anx_out <- matches[1,]
  } else {
    anx_out <- anx_out %>% add_row(matches[1,])
  }
}

anx_base <- select(anx_out,sjlid,anxdx)


# depression
dep <- behav[!is.na(behav$depdx),]

# select depression data closest to baseline neuro
test3 <- dep %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$datecomp)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_base <- grep(sid,base_neuro$sjlid, value = FALSE)
  neuro <- base_neuro[neuro_base,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    dep_out <- matches[1,]
  } else {
    dep_out <- dep_out %>% add_row(matches[1,])
  }
}

dep_base <- select(dep_out,sjlid,depdx)

merge1 <- merge(grd_base,anx_base, all=TRUE)
behav_base <- merge(merge1, dep_base, all=TRUE)


# behavior for recent neuro
# select grade data closest to baseline neuro
test3 <- grd %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$datecomp)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_recent <- grep(sid,recent_neuro$sjlid, value = FALSE)
  neuro <- recent_neuro[neuro_recent,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    grd_out <- matches[1,]
  } else {
    grd_out <- grd_out %>% add_row(matches[1,])
  }
}

grd_recent <- select(grd_out,sjlid,grade)


# select anxiety data closest to baseline neuro
test3 <- anx %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$datecomp)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_recent <- grep(sid,recent_neuro$sjlid, value = FALSE)
  neuro <- recent_neuro[neuro_recent,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    anx_out <- matches[1,]
  } else {
    anx_out <- anx_out %>% add_row(matches[1,])
  }
}

anx_recent <- select(anx_out,sjlid,anxdx)

# select depression data closest to baseline neuro
test3 <- dep %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$datecomp)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_recent <- grep(sid,recent_neuro$sjlid, value = FALSE)
  neuro <- recent_neuro[neuro_recent,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    dep_out <- matches[1,]
  } else {
    dep_out <- dep_out %>% add_row(matches[1,])
  }
}

dep_recent <- select(dep_out,sjlid,depdx)

merge1 <- merge(grd_recent,anx_recent, all=TRUE)
behav_recent <- merge(merge1, dep_recent, all=TRUE)






# physical activity
# adults
ahh <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Survey Data/adult_healthhabits.sas7bdat")
ahh2 <- select(ahh,SJLIFEID,agesurvey,datecomp,vpadays,vpamin,mpadays,mpamin)

# min of vigorous activity each week
ahh2$vminwk <- NULL
ahh2$vminwk <- ahh2$vpadays*ahh2$vpamin

# min of moderate activity each week
ahh2$mminwk <- NULL
ahh2$mminwk <- ahh2$mpadays*ahh2$mpamin

# total minutes of moderate and vig activity
ahh2$tminwk <- NULL
ahh2$tminwk <- rowSums(ahh2[,c("mminwk","vminwk")], na.rm=TRUE)

# physical activity: 1 = yes, 0 = no
ahh2$phyact <- ifelse(ahh2$vminwk >= 75 | ahh2$tminwk >= 150, 1, 0)
ahh2$phyact[is.na(ahh2$vminwk) & ahh2$tminwk >= 150] <- 1
ahh2$phyact[is.na(ahh2$vminwk) & ahh2$tminwk < 150 & ahh2$tminwk > 0] <- 0

# adolescent
khh <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Survey Data/adolescent_healthhabits.sas7bdat")
khh2 <- select(khh,SJLIFEID,agesurvey,datecomp,vpadays,vpamin,mpadays,mpamin)

# min of vigorous activity each week
khh2$vminwk <- NULL
khh2$vminwk <- khh2$vpadays*khh2$vpamin

# min of moderate activity each week
khh2$mminwk <- NULL
khh2$mminwk <- khh2$mpadays*khh2$mpamin

# total minutes of moderate and vig activity
khh2$tminwk <- NULL
khh2$tminwk <- rowSums(khh2[,c("mminwk","vminwk")], na.rm=TRUE)

# physical activity: 1 = yes, 0 = no
khh2$phyact <- ifelse(khh2$vminwk >= 75 | khh2$tminwk >= 150, 1, 0)
khh2$phyact[is.na(khh2$vminwk) & khh2$tminwk >= 150] <- 1
khh2$phyact[is.na(khh2$vminwk) & khh2$tminwk < 150 & khh2$tminwk > 0] <- 0

# abbreviated
abbrev <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Survey Data/abbreviated.sas7bdat")
abb2 <- select(abbrev,SJLIFEID,datecomp,mpamin,mpadays,vpadays,vpamin)

# min of vigorous activity each week
abb2$vminwk <- NULL
abb2$vminwk <- abb2$vpadays*abb2$vpamin

# min of moderate activity each week
abb2$mminwk <- NULL
abb2$mminwk <- abb2$mpadays*abb2$mpamin

# total minutes of moderate and vig activity
abb2$tminwk <- NULL
abb2$tminwk <- rowSums(abb2[,c("mminwk","vminwk")], na.rm=TRUE)

# physical activity: 1 = yes, 0 = no
abb2$phyact <- ifelse(abb2$vminwk >= 75 | abb2$tminwk >= 150, 1, 0)
abb2$phyact[is.na(abb2$vminwk) & abb2$tminwk >= 150] <- 1
abb2$phyact[is.na(abb2$vminwk) & abb2$tminwk < 150 & abb2$tminwk > 0] <- 0

abb2$agesurvey <- NA


merge1 <- rbind(ahh2,khh2)
phyact <- rbind(merge1,abb2)

# select
phyact <- phyact[!is.na(phyact$phyact),]
colnames(phyact) <- c("sjlid","agesurvey","datecomp","vpadays","vpamin","mpadays","mpamin","vminwk","mminwk","tminwk","phyact")

# select physical activity data closest to baseline neuro
test3 <- phyact %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$datecomp)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_base <- grep(sid,base_neuro$sjlid, value = FALSE)
  neuro <- base_neuro[neuro_base,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    pa_out <- matches[1,]
  } else {
    pa_out <- pa_out %>% add_row(matches[1,])
  }
}

pa_base <- select(pa_out,sjlid,phyact)


# recent
for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_recent <- grep(sid,recent_neuro$sjlid, value = FALSE)
  neuro <- recent_neuro[neuro_recent,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    pa_out <- matches[1,]
  } else {
    pa_out <- pa_out %>% add_row(matches[1,])
  }
}

pa_recent <- select(pa_out,sjlid,phyact)

# BMI 
func <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/function_combo_basic.sas7bdat")
func2 <- select(func,sjlid,DateVisitStart,BMIadj)
BMI <- func2[!is.na(func2$BMIadj),]

# select BMI data closest to baseline neuro
test3 <- BMI %>% group_by(sjlid) 
test3$datecomp <- as.Date(test3$DateVisitStart)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_base <- grep(sid,base_neuro$sjlid, value = FALSE)
  neuro <- base_neuro[neuro_base,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    bmi_out <- matches[1,]
  } else {
    bmi_out <- bmi_out %>% add_row(matches[1,])
  }
}

bmi_base <- select(bmi_out,sjlid,BMIadj)


# recent
for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_recent <- grep(sid,recent_neuro$sjlid, value = FALSE)
  neuro <- recent_neuro[neuro_recent,]
  matches$diff = abs(difftime(neuro$doe[1],matches$datecomp, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    bmi_out <- matches[1,]
  } else {
    bmi_out <- bmi_out %>% add_row(matches[1,])
  }
}

bmi_recent <- select(bmi_out,sjlid,BMIadj)


# pkv02 - need to request
func <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Neurocognition/functional.sas7bdat")

names(func)[names(func) == 'SJLID'] <- 'sjlid'

# remove all missing
func2 <- func[!is.na(func$PKVO2),] #| !is.na(func$SXMWHR_P) | !is.na(func$SXMWRPE_2) | !is.na(func$SXMWHR_2),]

# select closest to baseline
# select physical activity data closest to baseline neuro
test3 <- func2 %>% group_by(sjlid) 
test3$assmntdate <- as.Date(test3$assmntdate)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_base <- grep(sid,base_neuro$sjlid, value = FALSE)
  neuro <- base_neuro[neuro_base,]
  matches$diff = abs(difftime(neuro$doe[1],matches$assmntdate, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    func_out <- matches[1,]
  } else {
    func_out <- func_out %>% add_row(matches[1,])
  }
}

func_base <- select(func_out,sjlid,PKVO2)

# recent
for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_recent <- grep(sid,recent_neuro$sjlid, value = FALSE)
  neuro <- recent_neuro[neuro_recent,]
  matches$diff = abs(difftime(neuro$doe[1],matches$assmntdate, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    func_out <- matches[1,]
  } else {
    func_out <- func_out %>% add_row(matches[1,])
  }
}

func_recent <- select(func_out,sjlid,PKVO2)


# create length of follow-up variable for recent neurocognitive measures compared to base
recent_neuro$followup <- abs(difftime(recent_neuro$doe[1],base_neuro$doe, units = 'days'))

# merge data
merge1 <- merge(chemo2,base_lab, all = TRUE)
merge2 <- merge(merge1,rt2, all = TRUE)
merge3 <- merge(merge2,surg_out, all = TRUE)
merge4 <- merge(merge3,demo2, all = TRUE)
merge5 <- merge(merge4,diag3, all = TRUE)
merge6 <- merge(merge5,ctcae, all = TRUE)

base1 <- merge(merge6,smk_base, all = TRUE)
rec1 <- merge(merge6,smk_recent, all = TRUE)

base2 <- merge(base1,behav_base, all = TRUE)
rec2 <- merge(rec1,behav_recent, all = TRUE)

base3 <- merge(base2,pa_base, all = TRUE)
rec3 <- merge(rec2,pa_recent, all = TRUE)

base4 <- merge(base3,bmi_base, all = TRUE)
rec4 <- merge(rec3,bmi_recent, all = TRUE)

base5 <- merge(base4,func_base, all = TRUE)
rec5 <- merge(rec4,func_recent, all = TRUE)

baseline <- merge(base5,base_neuro, all.y = TRUE)
recent <- merge(rec5,recent_neuro, all.y = TRUE)


write.csv(baseline,"baseline_pheno_updated.csv", quote = FALSE, row.names = FALSE)
write.csv(recent,"recent_pheno_updated.csv", quote = FALSE, row.names = FALSE)


# get mrns
sjlid <- read.table("sjlids.txt", header = FALSE)
colnames(sjlid) <- 'sjlid'

demo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
demo <- select(demo,sjlid,mrn)

out <- merge(sjlid,demo)
out <- select(out,mrn)

write.csv(out,"mrns.csv", quote = FALSE, row.names = FALSE)



# read in new variable - CPTHitRTZs
new <- read_excel("addvariable.xlsx")
new2 <- merge(new,demo)
new2 <- select(new2,sjlid,CPTHitRTZs,doe)

new2 <- new2[!is.na(new2$CPTHitRTZs),]

# select smoking data closest to baseline neuro
test3 <- new2 %>% group_by(sjlid) 
test3$doe <- as.Date(test3$doe)

counts <- test3 %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


### matching base_neuro date
base_neuro <- read.csv("baseline_pheno_updated.csv", header = TRUE)
base_neuro <- select(base_neuro,sjlid,doe)

for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_base <- grep(sid,base_neuro$sjlid, value = FALSE)
  neuro <- base_neuro[neuro_base,]
  matches$diff = abs(difftime(neuro$doe[1],matches$doe, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    func_out <- matches[1,]
  } else {
    func_out <- func_out %>% add_row(matches[1,])
  }
}

func_out <- select(func_out,sjlid,CPTHitRTZs)

done <- merge(base_neuro, func_out, all.X = TRUE)

write.csv(done,"base_neuro_withCPTHitRTZs.csv", quote = FALSE, row.names = FALSE)


# udpate recent

recent_neuro <- read.csv("recent_pheno_updated.csv", header = TRUE)
recent_neuro <- select(recent_neuro,sjlid,doe)

for (i in 1:length(counts$n)) {
  # if counts = 1, only one measure 
  cat(i)
  sid <- counts$sjlid[i]
  filtered <- grep(sid,test3$sjlid, value=FALSE)
  matches <- test3[filtered,]
  matches <- matches %>% ungroup()
  neuro_recent <- grep(sid,recent_neuro$sjlid, value = FALSE)
  neuro <- recent_neuro[neuro_recent,]
  matches$diff = abs(difftime(neuro$doe[1],matches$doe, units = 'days'))
  matches <- matches %>% arrange(matches$diff)
  if ( i == 1) {
    func_out <- matches[1,]
  } else {
    func_out <- func_out %>% add_row(matches[1,])
  }
}

func_out <- select(func_out,sjlid,CPTHitRTZs)

done <- merge(recent_neuro,func_out, all.X=TRUE)
write.csv(done,"recent_neuro_withCPTHitRTZs.csv",quote = FALSE, row.names = FALSE)



base <- read.table("./assoc/base_pheno_withadmix.new.txt", header = TRUE)

colnames(func_out) <- c('FID','CPTHitRTZs')

update <- merge(base,func_out, all.x = TRUE)

write.table(update,"pheno_newest.txt", quote = FALSE, row.names = FALSE)


## add mrn
base <- read.csv("base_neuro_withCPTHitRTZs.csv", header = TRUE)

demo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
demo <- select(demo,sjlid,mrn)

out <- merge(demo,base)
write.csv(out,"base_neuro_withmrn.csv", quote = FALSE, row.names = FALSE)

recent <- read.csv("recent_neuro_withCPTHitRTZs.csv", header = TRUE)

out_rec <- merge(demo,recent)
write.csv(out_rec,"recent_neuro_withmrn.csv",quote = FALSE, row.names = FALSE)
