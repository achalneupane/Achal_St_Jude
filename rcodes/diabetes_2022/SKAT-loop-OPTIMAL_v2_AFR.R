#!/usr/bin/Rscript --vanilla --slave

## Achal-2022-09-14


### Script to run SKAT over a set of genes WITHIN A CHROMOSOME DIRECTORY
### USAGE: ./SKAT-unrelated-loop-OPTIMAL.R  ${WORKDIR} ${CHR} ${RAW} ${SET} ${COVARS}
library(SKAT)
library(tools)
library(data.table)
library(rrapply)
# Define arguments
# args <- commandArgs(TRUE)

workdir <- "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/gene-based-analysis/AFR/"
# workdir<-args[1]


covarsfile <- "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/SJLIFE_T2D_GWAS_AFR.pheno"
# covarsfile<-args[5]
# t2d,agedx,gender,age_last_visit,BMIadj,aa_class_dose_5,maxabdrtdose,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10
## LOAD COVARS file
cat(paste0("\nUploading covariates file ",covarsfile,"\n"))
covars <-read.table(covarsfile, head=T, check.names=F)

attach('Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA')
PHENO.ANY_SN <- PHENO.ANY_SN
detach('file:Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA')


sum(covars$IID %in% PHENO.ANY_SN$sjlid)
# 3113
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Could you please re-evaluate the top findings from this RV analysis with genes with SKAT-O P<1e-3 in treatment-stratified samples, per the following:; Re: 5
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1.a. Cutoff: for exposed, define as >200 cGy
# 1.b. Cutoff: same as above, but for RT cutoffs at >=500 cGy vs <500, >=1000 cGy vs <1000, >=1500 cGy vs
# 1.c Cutoff : same as above, but abdominal OR pelvic RT exposed
# abdominalRT dose
covars$maxabdrtdose <- PHENO.ANY_SN$maxabdrtdose[match(covars$IID, PHENO.ANY_SN$sjlid)]
# abdominalRT_YN 
covars$maxabdrtdose.exposed_more_than_200cGy_YN <- ifelse(covars$maxabdrtdose > 200, 1, 0)
covars$maxabdrtdose.exposed_500cGy_or_higher_YN <- ifelse(covars$maxabdrtdose >= 500, 1, 0)
covars$maxabdrtdose.exposed_1000cGy_or_higher_YN <- ifelse(covars$maxabdrtdose >= 1000, 1, 0)
covars$maxabdrtdose.exposed_1500cGy_or_higher_YN <- ifelse(covars$maxabdrtdose >= 1500, 1, 0)
covars$maxabdrtdose.exposed_2000cGy_or_higher_YN <- ifelse(covars$maxabdrtdose >= 2000, 1, 0)

# PelvisRT dose
covars$maxpelvisrtdose <- PHENO.ANY_SN$maxpelvisrtdose[match(covars$IID, PHENO.ANY_SN$sjlid)]
# PelvisRT_YN
covars$maxpelvisrtdose.exposed_more_than_200cGy_YN <- ifelse(covars$maxpelvisrtdose > 200, 1, 0)
covars$maxpelvisrtdose.exposed_500cGy_or_higher_YN <- ifelse(covars$maxpelvisrtdose >= 500, 1, 0)
covars$maxpelvisrtdose.exposed_1000cGy_or_higher_YN <- ifelse(covars$maxpelvisrtdose >= 1000, 1, 0)
covars$maxpelvisrtdose.exposed_1500cGy_or_higher_YN <- ifelse(covars$maxpelvisrtdose >= 1500, 1, 0)
covars$maxpelvisrtdose.exposed_2000cGy_or_higher_YN <- ifelse(covars$maxpelvisrtdose >= 2000, 1, 0)

# Abdominal or Pelvic RT exposed
covars$abd_OR_pelvis_exposed_more_than_200cGY_YN <- ifelse(covars$maxabdrtdose.exposed_more_than_200cGy_YN == 1 |covars$maxpelvisrtdose.exposed_more_than_200cGy_YN == 1, 1,0)
covars$abd_OR_pelvis.exposed_500cGy_or_higher_YN <- ifelse(covars$maxabdrtdose.exposed_500cGy_or_higher_YN == 1 |covars$maxpelvisrtdose.exposed_500cGy_or_higher_YN == 1, 1,0)
covars$abd_OR_pelvis.exposed_1000cGy_or_higher_YN <- ifelse(covars$maxabdrtdose.exposed_1000cGy_or_higher_YN == 1 |covars$maxpelvisrtdose.exposed_1000cGy_or_higher_YN == 1, 1,0)
covars$abd_OR_pelvis.exposed_1500cGy_or_higher_YN <- ifelse(covars$maxabdrtdose.exposed_1500cGy_or_higher_YN == 1 |covars$maxpelvisrtdose.exposed_1500cGy_or_higher_YN == 1, 1,0)
covars$abd_OR_pelvis.exposed_2000cGy_or_higher_YN <- ifelse(covars$maxabdrtdose.exposed_2000cGy_or_higher_YN == 1 |covars$maxpelvisrtdose.exposed_2000cGy_or_higher_YN == 1, 1,0)

# -------------------
# 2 alkylating agents
# -------------------
# AA within 5 years of primary cancer
covars$aa_class_dose_5 <- PHENO.ANY_SN$aa_class_dose_5[match(covars$IID, PHENO.ANY_SN$sjlid)]
covars$aa_class_dose_5_YN <- ifelse(covars$aa_class_dose_5 > 0, 1, 0 )
covars$aa_class_dose_5_4000_or_higher_YN <- ifelse(covars$aa_class_dose_5 >= 4000, 1, 0 )
# AA Any
covars$aa_class_dose_any <- PHENO.ANY_SN$aa_class_dose_any[match(covars$IID, PHENO.ANY_SN$sjlid)]
covars$aa_class_dose_any_YN <- ifelse(covars$aa_class_dose_any > 0, 1, 0 )
covars$aa_class_dose_any_4000_or_higher_YN <- ifelse(covars$aa_class_dose_any >= 4000, 1, 0 )

# --------------------------------------------------------------------------------------------
# 3. Could you please repeat step 5 but with cases defined using CTCAE-graded DM  >=3?; Re: 6
# --------------------------------------------------------------------------------------------
attach("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/CTCAE Abnormal Glucose Metabolism Rows With Transient Hyperglycemia-Mapping Rows Removed For Yadav's Population (n=4747, rows=9349).Rda")
diabetes29.ctcae <- diabetes29.ctcae
detach("file:Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/CTCAE Abnormal Glucose Metabolism Rows With Transient Hyperglycemia-Mapping Rows Removed For Yadav's Population (n=4747, rows=9349).Rda")
dim(diabetes29.ctcae)

## Extract samples with the most recent grade dates
diabetes29.ctcae <- diabetes29.ctcae[tapply(1:nrow(diabetes29.ctcae),diabetes29.ctcae$sjlid,function(ii) ii[which.max(diabetes29.ctcae$grade.date[ii])]),]
sum(covars$IID %in% covars$IID)
# 575

covars$ctcae_grade <- diabetes29.ctcae$grade[match(covars$IID, diabetes29.ctcae$sjlid)]
covars$ctcae_grade[covars$ctcae_grade == -9 ] <- NA
covars$ctcae_grad_3_or_higher_YN <- ifelse(covars$ctcae_grade >= 3, 1, 0)
# ------------------------------------------------------------------------------
# 4. Could you please use the treatment/case definition parameters set for 5 and
# 6 above but test for interaction (logit[y] ~ covariates + RV_burden +
# RV_burden * treatment), changing the definitions of "covariates" based on how
# you defined the treatment set in 5
# ------------------------------------------------------------------------------

################################################################
## 1. Running SKAT, BURDEN and SKAT-O using covars as in GWAS ##
################################################################

com.data <- covars
chrom <- c(1:17, 19:22)
for (j in 1:length(chrom)){
print(paste0("Doing chr", j))
chr<-chrom[j]
# chr<-args[2]

raw <- "AFR_diabetes_chrALL.dat-chr"
# raw<-args[3]

set<-"geneset-AFR"
# set<-args[4]

# set up directories:
chrdir<-paste0(workdir,"/chr", chr)
# rawfile<-paste0(chrdir,"/",raw,chr,".raw")
rawfile<-paste0(chrdir,"/",raw,chr,".raw")
setwd<-paste0(workdir, "/chr",chr)

#### LOAD RAW file
cat(paste0("\nUploading raw file ",rawfile," \n"))
system.time(raw <- read.table(rawfile, header=T, check.names=F) )
cat(dim(raw))

# arrange header of raw file
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER
raw <- raw[grep("chr|IID", colnames(raw))]
target_vec <- covars$IID
raw <- raw[match(target_vec, raw$IID),]
if(sum(com.data$IID != raw$IID) == 0){
com.data <- cbind.data.frame(com.data,raw[-1])
}
}


# #ARRANGE PHENO and GENDER so values are within  [0,1]
# #GENDER
# com.data$gender[com.data$gender==1]=0
# com.data$gender[com.data$gender==2]=1
# com.data$gender[com.data$gender==-9]=NA

#Pheno
com.data$t2d[com.data$t2d==-9]=NA
com.data$t2d[com.data$t2d==1]=0
com.data$t2d[com.data$t2d==2]=1


## Create gene list
geneLIST <- list()
geneLISTALL <- list()
extension<-".gene"
##  Create variant set by gene
datafiles<- Sys.glob(paste("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/gene-based-analysis/AFR/chr*/geneset-AFR","/*",extension,sep=""))
for (i in 1:length(datafiles)) {
  print(i)
  #print(datafiles[i])
  genename=basename(file_path_sans_ext(datafiles[i]))
  print(genename)
  
  ### GET GENO DATA
  SET<-read.table(datafiles[i],head=F, as.is=T, check.names=F,sep="\n")
  # geneLIST[i] <- SETlist 
  geneLIST[i] <- SET
  names(geneLIST)[i] <- genename
}

# rm(covars)

# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//diabetes/gene-based-analysis/rare_variant_analysis_AFR.Rdata")

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//diabetes/gene-based-analysis/rare_variant_analysis_AFR.Rdata")

#####################
## START LOOP PART ##
#####################

## Initialize values:
cat("\nInitialize values for table generation:\nGENE\tSNPID\tN0\tN1\tSKAT\tSKATO\tSKAT_125\n")
genelist<-"GENE"
SNPID<-"SNPID"
N0<-"N0"
N1<-"N1"
SKAT<-"SKAT"
SKATO<-"SKATO"
BURDEN<-"BURDEN"
### START loop over list of genes
for (i in 1:length(geneLIST)) {
  print(i)
  genename=names(geneLIST[i])
  print(genename)
  
  ### GET GENO DATA
  SETlist <- unname(rrapply(geneLIST[i],  how = 'unlist'))
  snpID<-paste(SETlist, collapse = ";")
  SETgeno<-as.matrix(com.data[,SETlist])
  
  ### Define y
  y=as.matrix(com.data$t2d)
  table(y)
  
  # Create RV_Burden variable
  com.data$RV_Burden <- rowSums(SETgeno)
  
  # MODEL adjusted by GWAS covars
  covariates = c("agedx","gender","age_last_visit","BMIadj","aa_class_dose_5","maxabdrtdose","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  index= which(colnames(com.data)%in%covariates)
  Xm2=as.matrix(com.data[,index])
  obj<-SKAT_Null_Model(y ~ Xm2,out_type="D")
  n0<-SKAT(SETgeno, obj, kernel = "linear.weighted")$param$n.marker	# n input SNPs
  n1<-SKAT(SETgeno, obj, kernel = "linear.weighted")$param$n.marker.test # n ouptut SNPs
  skat<-SKAT(SETgeno, obj, kernel = "linear.weighted")$p.value
  skato<-SKAT(SETgeno, obj, kernel = "linear.weighted", method="optimal.adj")$p.value
  ## Burden can be stated in two different ways
  # burden<-SKAT(SETgeno, obj, kernel = "linear.weighted", method="Burden")$p.value
  burden<-SKAT(SETgeno, obj, kernel = "linear.weighted", method="davies", r.corr =1)$p.value
  
  ## add values to lists
  genelist<-rbind(genelist,genename)
  SNPID<-rbind(SNPID,snpID)
  N0<-rbind(N0,n0)
  N1<-rbind(N1,n1)
  
  SKAT<-rbind(SKAT,skat)
  SKATO<-rbind(SKATO,skato)
  BURDEN <- rbind(BURDEN, burden)
}


chrALL <- cbind.data.frame(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
colnames(chrALL) <- chrALL[1,]
chrALL <- chrALL[-1,]
chrALL[,c("SKAT", "SKATO", "BURDEN")] <- sapply(chrALL[,c("SKAT", "SKATO", "BURDEN")], as.numeric)

# sum(chrALL$SKATO < 1e-3)
# 

write.table(chrALL, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-AFR_chrALL_SKAT.pval", sep="\t", col.names=T, row.names=F, quote=F)


##################################################
## Re-evaluate the top finding (question no. 5) ##
##################################################
# Could you please re-evaluate the top findings from this RV analysis with genes with SKAT-O P<1e-3 in treatment-stratified samples, per the following:
# abdominal RT exposed (for exposed, define as >200 cGy)
# pelvic RT exposed
# abdominal OR pelvic RT exposed
# same as above, but for RT cutoffs at >=500 cGy vs <500, >=1000 cGy vs <1000, >=1500 cGy vs <1500, >=2000 cGy vs <2000
# alkylating agents exposed (none vs any)
# alkylating agents, using dose cutoffs at >=4000 mg/m2 vs <4000 mg/m2

genelist<-"GENE"
SNPID<-"SNPID"
N0<-"N0" # Number of input SNPs
N1<-"N1" # Number of tested SNPs
SKAT<-"SKAT"
SKATO<-"SKATO"
BURDEN<-"BURDEN"


wanted.genes <- chrALL$GENE[chrALL$SKATO < 1e-3]
geneLIST.wanted <- geneLIST[names(geneLIST) %in% wanted.genes]



# a. abdominal RT exposed (for exposed, define as >200 cGy)
### START loop over list of genes
for (i in 1:length(geneLIST.wanted)) {
  print(i)
  genename=names(geneLIST.wanted[i])
  print(genename)
  
  ### GET GENO DATA
  SETlist <- unname(rrapply(geneLIST.wanted[i],  how = 'unlist'))
  snpID<-paste(SETlist, collapse = ";")
  SETgeno<-as.matrix(com.data[,SETlist])
  
  ### Define y
  y=as.matrix(com.data$t2d)
  table(y)
  
  # Create RV_Burden variable
  com.data$RV_Burden <- rowSums(SETgeno)
  
  # MODEL adjusted
  covariates = c("agedx","gender","age_last_visit","BMIadj","aa_class_dose_5","maxabdrtdose","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

  index= which(colnames(com.data)%in%covariates)
  Xm2=as.matrix(com.data[,index])
  obj<-SKAT_Null_Model(y ~ Xm2,out_type="D")
  n0<-SKAT(SETgeno, obj, kernel = "linear.weighted")$param$n.marker	# n input SNPs
  n1<-SKAT(SETgeno, obj, kernel = "linear.weighted")$param$n.marker.test # n ouptut SNPs
  skat<-SKAT(SETgeno, obj, kernel = "linear.weighted")$p.value
  skato<-SKAT(SETgeno, obj, kernel = "linear.weighted", method="optimal.adj")$p.value
  ## Burden can be stated in two different ways
  # burden<-SKAT(SETgeno, obj, kernel = "linear.weighted", method="Burden")$p.value
  burden<-SKAT(SETgeno, obj, kernel = "linear.weighted", method="davies", r.corr =1)$p.value
  
  ## add values to lists
  genelist<-rbind(genelist,genename)
  SNPID<-rbind(SNPID,snpID)
  N0<-rbind(N0,n0)
  N1<-rbind(N1,n1)
  
  SKAT<-rbind(SKAT,skat)
  SKATO<-rbind(SKATO,skato)
  BURDEN <- rbind(BURDEN, burden)
}


chrALL <- cbind.data.frame(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
colnames(chrALL) <- chrALL[1,]
chrALL <- chrALL[-1,]
chrALL[,c("SKAT", "SKATO", "BURDEN")] <- sapply(chrALL[,c("SKAT", "SKATO", "BURDEN")], as.numeric)









# Calculate RV_Burden
com.data$RV_Burden <- rowSums(SETgeno)
# print(max(com.data$RV_Burden, na.rm = T))
