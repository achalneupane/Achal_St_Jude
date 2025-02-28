### pull variables for HTN GWAS

library(haven)
library(dplyr)

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/hypertension/")

# read in IDs
idsEA <- read.table("keep_SJLIFE_EA.txt", header = FALSE) 
colnames(idsEA) <- c('sjlid','remove')

idsEA <- select(idsEA,sjlid)

idsAA <- read.table("keep_SJLIFE_AA.txt", header = FALSE)
colnames(idsAA) <- c('sjlid','remove')

idsAA <- select(idsAA,sjlid)

all <- read.table("ids_all.txt", header = FALSE)
colnames(all) <- 'sjlid'

### demographics
demo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
demo2 <- select(demo,sjlid,gender)

demoAA <- merge(idsAA,demo2)
demoEA <- merge(idsEA,demo2)
demoall <- merge(all,demo2)

# diagnosis
diag <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat")
diag <- diag[diag$primdx == 1,]
diag2 <- select(diag,sjlid,agedx)

diagAA <- merge(idsAA,diag2)
diagEA <- merge(idsEA,diag2)
diagall <- merge(all,diag2)

# last contact
lst <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Tracking Data/lstcondt.sas7bdat")
lst2 <- select(lst,sjlid,agelstcontact)

lstAA <- merge(idsAA,lst2)
lstEA <- merge(idsEA,lst2)
lstall <- merge(all,lst2)



# radiation # abdomen and HPA (seg2?) ################ check with yadav
rad <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat")
rad2 <- select(rad,sjlid,maxabdrtdose,TBIDose)

radAA <- merge(idsAA,rad2)
radEA <- merge(idsEA,rad2)
radall <- merge(all,rad2)

# chemo 
chemo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat")
chemo2 <- select(chemo,sjlid,alkylating_dose_5,platinum_ccss_dose_5)

chemoAA <- merge(idsAA,chemo2)
chemoEA <- merge(idsEA,chemo2)
chemoall <- merge(all,chemo2)

# BP
milli <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/milli_measures.sas7bdat")
milli2 <- select(milli,sjlid,labdt,DBP_MAX,DBP_MIN,SBP_MAX,SBP_MIN)

milli2 <- milli2[milli2$DBP_MAX >= 1,]

# select baseline
# sort by date - pick baseline for each individual
test <- milli2 %>% group_by(sjlid) 

counts <- test %>% summarise(n=n())

counts$n <- as.numeric(counts$n)

col_names <- c("sjlid","labdt","DBP_MAX","DBP_MIN","SBP_MAX","SBP_MIN")
base_bp <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
colnames(base_bp) <- col_names
base_bp$sjlid <- as.character(base_bp$sjlid)
base_bp$labdt <- as.Date(base_bp$labdt)
base_bp$DBP_MAX <- as.character(base_bp$DBP_MAX)
base_bp$DBP_MIN <- as.character(base_bp$DBP_MIN)
base_bp$SBP_MAX <- as.character(base_bp$SBP_MAX)
base_bp$SBP_MIN <- as.character(base_bp$SBP_MIN)



for (i in 1:length(counts$n)) {
  cat(i)
  sid <- counts$sjlid[i]
  # test3 <- test3 %>% ungroup()
  filtered <- grep(sid,test$sjlid, value=FALSE)
  matches <- test[filtered,]
  matches <- matches %>% arrange(labdt)
  
  base_bp <- base_bp %>% add_row(matches[1,])
}


# recent bp measures

col_names <- c("sjlid","labdt","DBP_MAX","DBP_MIN","SBP_MAX","SBP_MIN")
recent_bp <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
colnames(recent_bp) <- col_names
recent_bp$sjlid <- as.character(recent_bp$sjlid)
recent_bp$labdt <- as.Date(recent_bp$labdt)
recent_bp$DBP_MAX <- as.character(recent_bp$DBP_MAX)
recent_bp$DBP_MIN <- as.character(recent_bp$DBP_MIN)
recent_bp$SBP_MAX <- as.character(recent_bp$SBP_MAX)
recent_bp$SBP_MIN <- as.character(recent_bp$SBP_MIN)



for (i in 1:length(counts$n)) {
  cat(i)
  sid <- counts$sjlid[i]
  # test3 <- test3 %>% ungroup()
  filtered <- grep(sid,test$sjlid, value=FALSE)
  matches <- test[filtered,]
  matches <- matches %>% arrange(desc(labdt))
  
  recent_bp <- recent_bp %>% add_row(matches[1,])
}



base_bp <- base_bp[,-2]

base_bp$DBP_MAX <- as.numeric(base_bp$DBP_MAX)
base_bp$DBP_MIN <- as.numeric(base_bp$DBP_MIN)
base_bp$SBP_MAX <- as.numeric(base_bp$SBP_MAX)
base_bp$SBP_MIN <- as.numeric(base_bp$SBP_MIN)

base_bp$PP_MIN <- (base_bp$SBP_MIN - base_bp$DBP_MIN)
base_bp$PP_MAX <- (base_bp$SBP_MAX - base_bp$DBP_MAX)


recent_bp <- recent_bp[,-2]

recent_bp$DBP_MAX <- as.numeric(recent_bp$DBP_MAX)
recent_bp$DBP_MIN <- as.numeric(recent_bp$DBP_MIN)
recent_bp$SBP_MAX <- as.numeric(recent_bp$SBP_MAX)
recent_bp$SBP_MIN <- as.numeric(recent_bp$SBP_MIN)

recent_bp$PP_MIN <- (recent_bp$SBP_MIN - recent_bp$DBP_MIN)
recent_bp$PP_MAX <- (recent_bp$SBP_MAX - recent_bp$DBP_MAX)

bpEA_base <- merge(idsEA,base_bp)
colnames(bpEA_base) <- c('sjlid','DBP_MAX_B','DBP_MIN_B','SBP_MAX_B','SBP_MIN_B',"PP_MIN_B","PP_MAX_B")
bpEA_rec <- merge(idsEA,recent_bp)
colnames(bpEA_rec) <- c('sjlid','DBP_MAX_R','DBP_MIN_R','SBP_MAX_R','SBP_MIN_R',"PP_MIN_R","PP_MAX_R")
bpEA <- merge(bpEA_base,bpEA_rec)

bpAA_base <- merge(idsAA,base_bp)
colnames(bpAA_base) <- c('sjlid','DBP_MAX_B','DBP_MIN_B','SBP_MAX_B','SBP_MIN_B',"PP_MIN_B","PP_MAX_B")
bpAA_rec <- merge(idsAA,recent_bp)
colnames(bpAA_rec) <- c('sjlid','DBP_MAX_R','DBP_MIN_R','SBP_MAX_R','SBP_MIN_R',"PP_MIN_R","PP_MAX_R")
bpAA <- merge(bpAA_base,bpAA_rec)


bpall_base <- merge(all,base_bp)
colnames(bpall_base) <- c('sjlid','DBP_MAX_B','DBP_MIN_B','SBP_MAX_B','SBP_MIN_B',"PP_MIN_B","PP_MAX_B")
bpall_rec <- merge(all,recent_bp)
colnames(bpall_rec) <- c('sjlid','DBP_MAX_R','DBP_MIN_R','SBP_MAX_R','SBP_MIN_R',"PP_MIN_R","PP_MAX_R")
bpall <- merge(bpall_base,bpall_rec)

### CTCAE
ctcae <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
hpt <- ctcae[ctcae$condition == "Hypertension",]
hpt2 <- select(hpt,sjlid,gradedt,grade)

# sort by grade - group by sjlid
test <- hpt2 %>% group_by(sjlid) 

counts <- test %>% summarise(n=n())

counts$n <- as.numeric(counts$n)


col_names <- c("sjlid","gradedt","grade")
grade <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
colnames(grade) <- col_names
grade$sjlid <- as.character(grade$sjlid)
grade$gradedt <- as.Date(grade$gradedt)
grade$grade <- as.numeric(grade$grade)



for (i in 1:length(counts$n)) {
  cat(i)
  sid <- counts$sjlid[i]
  # test3 <- test3 %>% ungroup()
  filtered <- grep(sid,test$sjlid, value=FALSE)
  matches <- test[filtered,]
  matches <- matches %>% arrange(desc(grade))
  
  grade <- grade %>% add_row(matches[1,])
}

grade$hpt <- ifelse(grade$grade >= 2,2,1)

grade <- select(grade,sjlid,hpt)

gradeEA <- merge(idsEA,grade)
gradeAA <- merge(idsAA,grade)
gradeall <- merge(all,grade)


EA1 <- merge(demoEA,diagEA)
EA2 <- merge(EA1,lstEA)
EA3 <- merge(EA2,radEA)
EA4 <- merge(EA3,chemoEA)
EA5 <- merge(EA4,bpEA, all.x = TRUE)
EA6 <- merge(EA5,gradeEA)

write.table(EA6,"pheno_EA_hp.txt", quote = FALSE, row.names = FALSE)


AA1 <- merge(demoAA,diagAA)
AA2 <- merge(AA1,lstAA)
AA3 <- merge(AA2,radAA)
AA4 <- merge(AA3,chemoAA)
AA5 <- merge(AA4,bpAA, all.x = TRUE)
AA6 <- merge(AA5,gradeAA)

write.table(AA6,"pheno_AA_hp.txt", quote = FALSE, row.names = FALSE)

all1 <- merge(demoall,diagall)
all2 <- merge(all1,lstall)
all3 <- merge(all2,radall)
all4 <- merge(all3,chemoall)
all5 <- merge(all4,bpall, all.x = TRUE)
all6 <- merge(all5,gradeall)

admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
colnames(admixture) <- c("sjlid","EUR","EAS","AFR")
admix <- select(admixture,sjlid,EAS,AFR)

out <- merge(all6,admix)

write.table(out,"pheno_all_hp.txt", quote = FALSE, row.names = FALSE)


AA <- read.table("pheno_AA_hp.txt", header = TRUE)

AAPCs <- read.table("AFR.eigenvec", header = TRUE)
AAPCs <- select(AAPCs,FID,PC1,PC2,PC3)
colnames(AAPCs) <- c('sjlid','PC1','PC2','PC3')

AAout <- merge(AA,AAPCs)

write.table(AAout,"pheno_AA_hp.txt", quote = FALSE, row.names = FALSE)


EA <- read.table("pheno_EA_hp.txt", header = TRUE)

EAPCs <- read.table("EUR.eigenvec", header = TRUE)
EAPCs <- select(EAPCs,FID,PC1,PC2,PC3)
colnames(EAPCs) <- c('sjlid','PC1','PC2','PC3')

EAout <- merge(EA,EAPCs)
write.table(EAout,"pheno_EA_hp.txt", quote = FALSE, row.names = FALSE)

# readd hpt
AFR <- read.table("pheno_AFR_hp.txt", header = TRUE)
AFR$hpt[AFR$grade >= 2] <- 2
AFR$hpt[AFR$grade == 0] <- 1
AFR$hpt[AFR$grade == 1] <- NA

write.table(AFR,"pheno_AFR_hp.txt", row.names = FALSE, quote = FALSE)

EUR <- read.table("pheno_EUR_hp.txt", header = TRUE)
EUR$hpt[EUR$grade >= 2] <- 2
EUR$hpt[EUR$grade == 0] <- 1
EUR$hpt[EUR$grade == 1] <- NA
write.table(EUR,"pheno_EUR_hp.txt", quote = FALSE, row.names = FALSE)

all <- read.table("pheno_all_hp.txt", header = TRUE)
all$hpt[all$grade >= 2] <- 2
all$hpt[all$grade == 0] <- 1
all$hpt[all$grade == 1] <- NA
write.table(all,"pheno_all_hp.txt", quote = FALSE, row.names = FALSE)




#### pheno --- add new pcs
EUR <- read.table("pheno_EUR_hp.txt", header = TRUE)
EUR <- EUR[,-(23:25)]
Epcs <- read.table("EUR.new.eigenvec", header = TRUE)
EPCS <- select(Epcs,FID,PC1,PC2,PC3)

EURout <- merge(EUR,EPCS)
write.table(EURout,"Pheno_EUR_hp_new.txt", quote = FALSE, row.names = FALSE)

AFR <- read.table("pheno_AFR_hp.txt", header = TRUE)
AFR <- AFR[,-(23:25)]

Apcs <- read.table("AFR.new.eigenvec", header = TRUE)
APCS <- select(Apcs,FID,PC1,PC2,PC3)
AFRout <- merge(AFR,APCS)
write.table(AFRout,"pheno_AFR_hp_new.txt", quote = FALSE, row.names = FALSE)


### adjust bp for medication - add 15 to SBP and 10 to DBP
EUR <- read.table("pheno_EUR_hp.txt", header = TRUE)

AFR <- read.table("pheno_AFR_hp.txt", header = TRUE)

all <- read.table("pheno_all_hp.txt", header = TRUE)

EUR$SBP_new <- ifelse(EUR$grade >=2, EUR$SBP_MAX_B + 15,EUR$SBP_MAX_B)
EUR$DBP_new <- ifelse(EUR$grade >=2, EUR$DBP_MAX_B + 10,EUR$DBP_MAX_B)
EUR$PP_new <- EUR$SBP_new - EUR$DBP_new


AFR$SBP_new <- ifelse(AFR$grade >=2, AFR$SBP_MAX_B + 15,AFR$SBP_MAX_B)
AFR$DBP_new <- ifelse(AFR$grade >=2, AFR$DBP_MAX_B + 10,AFR$DBP_MAX_B)
AFR$PP_new <- AFR$SBP_new - AFR$DBP_new

all$SBP_new <- ifelse(all$grade >=2, all$SBP_MAX_B + 15,all$SBP_MAX_B)
all$DBP_new <- ifelse(all$grade >=2, all$DBP_MAX_B + 10,all$DBP_MAX_B)
all$PP_new <- all$SBP_new - all$DBP_new

write.table(EUR, "pheno_EUR_fixed.txt", quote = FALSE, row.names = FALSE)
write.table(AFR, "pheno_AFR_fixed.txt", quote = FALSE, row.names = FALSE)
write.table(all,"pheno_all_fixed.txt", quote = FALSE, row.names = FALSE)
