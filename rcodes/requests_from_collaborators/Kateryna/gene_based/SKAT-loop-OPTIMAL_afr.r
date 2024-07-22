#!/usr/bin/Rscript --vanilla --slave

## Achal-2022-09-14


### Script to run SKAT over a set of genes WITHIN A CHROMOSOME DIRECTORY
### USAGE: ./SKAT-unrelated-loop-OPTIMAL.R  ${WORKDIR} ${CHR} ${RAW} ${SET} ${COVARS}
library(SKAT)
library(tools)
library(data.table)
# Define arguments
# args <- commandArgs(TRUE)

workdir <- "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/AFR/"
# workdir<-args[1]

chrALL <- {}
chrom <- 1:22
for (j in 1:length(chrom)){
print(paste0("Doing chr", j))
chr<-chrom[j]
# chr<-args[2]

raw <- "AFR_final-chr"
# raw<-args[3]

set<-"geneset-AFR"
# set<-args[4]

covarsfile <- "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/pheno/sjlife_afr_dox_only_pcs.pheno"
# covarsfile<-args[5]
# t2d,agedx,gender,age_last_visit,BMIadj,aa_class_dose_5,maxabdrtdose,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

# set up directories:
chrdir<-paste0(workdir,"/chr", chr)
# rawfile<-paste0(chrdir,"/",raw,chr,".raw")
rawfile<-paste0(chrdir,"/",raw,chr,".raw")
setwd<-paste0(workdir, "/chr",chr)

#### LOAD RAW file
cat(paste0("\nUploading raw file ",rawfile," \n"))
# system.time(raw <- read.table(rawfile, header=T, check.names=F) )

raw <- tryCatch({
  read.table(rawfile, header = TRUE, check.names = FALSE)
}, error = function(e) {
  message(paste("Skipping chromosome", chr, "due to error: File not found"))
  return(NULL)
})

if (is.null(raw)) {
  next
}

cat(dim(raw))




# arrange header of raw file
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

## LOAD COVARS file
cat(paste0("\nUploading covariates file ",covarsfile,"\n"))
covars <-read.table(covarsfile, head=T, check.names=F)

#MERGE1=merge(covars,raw, by="IID" ) #use this way and follwoing STEPS IF adding further covariants
com.data=merge(covars,raw, by="IID" )
cat("\n\ndimensions of com.data matrix are:\n")
dim(com.data)

# #ARRANGE PHENO and GENDER so values are within  [0,1]
# #GENDER
# com.data$gender[com.data$gender==1]=0
# com.data$gender[com.data$gender==2]=1

#Pheno
com.data$CMP2plus[com.data$CMP2plus==1]=0
com.data$CMP2plus[com.data$CMP2plus==2]=1

#####################
## START LOOP PART ##
#####################
## CHANGE DIR
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
  y=as.matrix(com.data$CMP2plus)
  table(y)
  
  # MODEL adjusted by SEX and fisrt three PCs
  covariates = c("agedx","gender","agelstcontact","anthra_jco_dose_any","Chest","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
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


colnames(chrALL) <- c("genelist","SNPID","N0","N1","SKAT","SKATO","BURDEN")

chrALL$SKAT_FDR <- p.adjust(as.numeric(chrALL$SKAT), method = "fdr")
chrALL$SKATO_FDR <- p.adjust(as.numeric(chrALL$SKATO), method = "fdr")
chrALL$BURDEN_FDR <- p.adjust(as.numeric(chrALL$BURDEN), method = "fdr")

write.table(chrALL, paste(workdir, set, "_chrALL","_SKAT_dox.pval",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
