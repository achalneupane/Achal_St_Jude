#!/usr/bin/Rscript --vanilla --slave

## Achal-2022-09-14


### Script to run SKAT over a set of genes WITHIN A CHROMOSOME DIRECTORY
### USAGE: ./SKAT-unrelated-loop-OPTIMAL.R  ${WORKDIR} ${CHR} ${RAW} ${SET} ${COVARS}
library(SKAT)
library(tools)
library(data.table)
# Define arguments
# args <- commandArgs(TRUE)

workdir <- "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/gene-based-analysis/EUR/"
# workdir<-args[1]


covarsfile <- "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/SJLIFE_T2D_GWAS_EUR.pheno"
# covarsfile<-args[5]
# t2d,agedx,gender,age_last_visit,BMIadj,aa_class_dose_5,maxabdrtdose,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10
## LOAD COVARS file
cat(paste0("\nUploading covariates file ",covarsfile,"\n"))
covars <-read.table(covarsfile, head=T, check.names=F)

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v2.RDATA")

sum(covars$IID %in% PHENO.ANY_SN$sjlid)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Could you please re-evaluate the top findings from this RV analysis with genes with SKAT-O P<1e-3 in treatment-stratified samples, per the following:; Re: 5
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1.a. Cutoff: for exposed, define as >200 cGy
# 1.b. Cutoff: same as above, but for RT cutoffs at >=500 cGy vs <500, >=1000 cGy vs <1000, >=1500 cGy vs

# abdominalRT dose
covars$maxabdrtdose <- PHENO.ANY_SN$maxabdrtdose[match(covars$IID, PHENO.ANY_SN$sjlid)]
# PelvisRT dose
covars$maxpelvisrtdose <- PHENO.ANY_SN$maxpelvisrtdose[match(covars$IID, PHENO.ANY_SN$sjlid)]

# -------------------
# 2 alkylating agents
# -------------------
# AA within 5 years of primary cancer
covars$aa_class_dose_5 <- PHENO.ANY_SN$aa_class_dose_5[match(covars$IID, PHENO.ANY_SN$sjlid)]
# AA Any
covars$aa_class_dose_any <- PHENO.ANY_SN$aa_class_dose_any[match(covars$IID, PHENO.ANY_SN$sjlid)]

# --------------------------------------------------------------------------------------------
# 3. Could you please repeat step 5 but with cases defined using CTCAE-graded DM  >=3?; Re: 6
# --------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 4. Could you please use the treatment/case definition parameters set for 5 and
# 6 above but test for interaction (logit[y] ~ covariates + RV_burden +
# RV_burden * treatment), changing the definitions of "covariates" based on how
# you defined the treatment set in 5
# ------------------------------------------------------------------------------


chrALL <- {}
chrom <- 1:22
for (j in 1:length(chrom)){
print(paste0("Doing chr", j))
chr<-chrom[j]
# chr<-args[2]

raw <- "EUR_diabetes_chrALL.dat-chr"
# raw<-args[3]

set<-"geneset-EUR"
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

#MERGE1=merge(covars,raw, by="IID" ) #use this way and follwoing STEPS IF adding further covariants
com.data=merge(covars,raw, by="IID" )
cat("\n\ndimensions of com.data matrix are:\n")
dim(com.data)

# #ARRANGE PHENO and GENDER so values are within  [0,1]
# #GENDER
# com.data$gender[com.data$gender==1]=0
# com.data$gender[com.data$gender==2]=1
# com.data$gender[com.data$gender==-9]=NA

#Pheno
com.data$t2d[com.data$t2d==-9]=NA
com.data$t2d[com.data$t2d==1]=0
com.data$t2d[com.data$t2d==2]=1

#####################
## START LOOP PART ##
#####################
## CHANGE DIR
#set<-"geneset-MAF-0.005"
dir<-paste0(chrdir,"/",set)
extension<-".gene"

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
cat(paste0("\nStart loop over list of genes in chr ",chr,"in geneset driectory",dir,"\n"))
datafiles<-Sys.glob(paste(dir,"/*",extension,sep=""))
for (i in 1:length(datafiles)) {
  print(i)
  #print(datafiles[i])
  genename=basename(file_path_sans_ext(datafiles[i]))
  print(genename)
  
  ### GET GENO DATA
  SET<-read.table(datafiles[i],head=F, as.is=T, check.names=F,sep="\n")
  SETlist<-c(SET$V1)
  # snpID<-SET[1,1]
  snpID<-paste(SETlist, collapse = ";")
  SETgeno<-as.matrix(com.data[,SETlist])
  
  ### Define y
  y=as.matrix(com.data$t2d)
  table(y)
  
  # MODEL adjusted by SEX and fisrt three PCs
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


chr.tmp <- cbind(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
colnames(chr.tmp) <- chr.tmp[1,]
chr.tmp <- chr.tmp[-1,]
chrALL <- rbind.data.frame(chrALL, chr.tmp)
}

write.table(chrALL, paste(workdir, set, "_chrALL","_SKAT.pval",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
