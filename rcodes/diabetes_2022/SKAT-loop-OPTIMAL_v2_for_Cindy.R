################################
## Create variables for Cindy ##
################################
# # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# # 1. Could you please re-evaluate the top findings from this RV analysis with genes with SKAT-O P<1e-3 in treatment-stratified samples, per the following:; Re: 5
# # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# # 1.a. Cutoff: for exposed, define as >200 cGy
# # 1.b. Cutoff: same as above, but for RT cutoffs at >=500 cGy vs <500, >=1000 cGy vs <1000, >=1500 cGy vs
# # 1.c Cutoff : same as above, but abdominal OR pelvic RT exposed
# # abdominalRT dose
# # abdominalRT_YN 
# covars$maxabdrtdose.exposed_more_than_200cGy_YN <- ifelse(covars$maxabdrtdose > 200, "Y", "N")
# covars$maxabdrtdose.exposed_500cGy_or_higher_YN <- ifelse(covars$maxabdrtdose >= 500, "Y", "N")
# covars$maxabdrtdose.exposed_1000cGy_or_higher_YN <- ifelse(covars$maxabdrtdose >= 1000, "Y", "N")
# covars$maxabdrtdose.exposed_1500cGy_or_higher_YN <- ifelse(covars$maxabdrtdose >= 1500, "Y", "N")
# covars$maxabdrtdose.exposed_2000cGy_or_higher_YN <- ifelse(covars$maxabdrtdose >= 2000, "Y", "N")
# 
# # PelvisRT dose
# covars$maxpelvisrtdose <- PHENO.ANY_SN$maxpelvisrtdose[match(covars$IID, PHENO.ANY_SN$sjlid)]
# # PelvisRT_YN
# covars$maxpelvisrtdose.exposed_more_than_200cGy_YN <- ifelse(covars$maxpelvisrtdose > 200, "Y", "N")
# covars$maxpelvisrtdose.exposed_500cGy_or_higher_YN <- ifelse(covars$maxpelvisrtdose >= 500, "Y", "N")
# covars$maxpelvisrtdose.exposed_1000cGy_or_higher_YN <- ifelse(covars$maxpelvisrtdose >= 1000, "Y", "N")
# covars$maxpelvisrtdose.exposed_1500cGy_or_higher_YN <- ifelse(covars$maxpelvisrtdose >= 1500, "Y", "N")
# covars$maxpelvisrtdose.exposed_2000cGy_or_higher_YN <- ifelse(covars$maxpelvisrtdose >= 2000, "Y", "N")
# 
# # Abdominal or Pelvic RT exposed
# covars$abd_OR_pelvis_exposed_more_than_200cGY_YN <- ifelse(covars$maxabdrtdose.exposed_more_than_200cGy_YN == "Y" |covars$maxpelvisrtdose.exposed_more_than_200cGy_YN == "Y", "Y","N")
# covars$abd_OR_pelvis.exposed_500cGy_or_higher_YN <- ifelse(covars$maxabdrtdose.exposed_500cGy_or_higher_YN == "Y" |covars$maxpelvisrtdose.exposed_500cGy_or_higher_YN == "Y", "Y","N")
# covars$abd_OR_pelvis.exposed_1000cGy_or_higher_YN <- ifelse(covars$maxabdrtdose.exposed_1000cGy_or_higher_YN == "Y" |covars$maxpelvisrtdose.exposed_1000cGy_or_higher_YN == "Y", "Y","N")
# covars$abd_OR_pelvis.exposed_1500cGy_or_higher_YN <- ifelse(covars$maxabdrtdose.exposed_1500cGy_or_higher_YN == "Y" |covars$maxpelvisrtdose.exposed_1500cGy_or_higher_YN == "Y", "Y","N")
# covars$abd_OR_pelvis.exposed_2000cGy_or_higher_YN <- ifelse(covars$maxabdrtdose.exposed_2000cGy_or_higher_YN == "Y" |covars$maxpelvisrtdose.exposed_2000cGy_or_higher_YN == "Y", "Y","N")
# 
# # -------------------
# # 2 alkylating agents
# # -------------------
# # AA within 5 years of primary cancer
# covars$aa_class_dose_5 <- PHENO.ANY_SN$aa_class_dose_5[match(covars$IID, PHENO.ANY_SN$sjlid)]
# covars$aa_class_dose_5_YN <- ifelse(covars$aa_class_dose_5 > 0, "Y", "N" )
# covars$aa_class_dose_5_4000_or_higher_YN <- ifelse(covars$aa_class_dose_5 >= 4000, "Y", "N" )
# # AA Any
# covars$aa_class_dose_any <- PHENO.ANY_SN$aa_class_dose_any[match(covars$IID, PHENO.ANY_SN$sjlid)]
# covars$aa_class_dose_any_YN <- ifelse(covars$aa_class_dose_any > 0, "Y", "N" )
# covars$aa_class_dose_any_4000_or_higher_YN <- ifelse(covars$aa_class_dose_any >= 4000, "Y", "N" )
# 
# # --------------------------------------------------------------------------------------------
# # 3. Could you please repeat step 5 but with cases defined using CTCAE-graded DM  >=3?; Re: 6
# # --------------------------------------------------------------------------------------------
# attach("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/CTCAE Abnormal Glucose Metabolism Rows With Transient Hyperglycemia-Mapping Rows Removed For Yadav's Population (n=4747, rows=9349).Rda")
# diabetes29.ctcae <- diabetes29.ctcae
# detach("file:Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/CTCAE Abnormal Glucose Metabolism Rows With Transient Hyperglycemia-Mapping Rows Removed For Yadav's Population (n=4747, rows=9349).Rda")
# dim(diabetes29.ctcae)
# 
# ## Extract samples with the most recent grade dates
# diabetes29.ctcae <- diabetes29.ctcae[tapply(1:nrow(diabetes29.ctcae),diabetes29.ctcae$sjlid,function(ii) ii[which.max(diabetes29.ctcae$grade.date[ii])]),]
# sum(covars$IID %in% covars$IID)
# # 3113
# 
# covars$ctcae_grade <- diabetes29.ctcae$grade[match(covars$IID, diabetes29.ctcae$sjlid)]
# covars$ctcae_grade[covars$ctcae_grade == -9 ] <- NA
# covars$ctcae_grad_3_or_higher_YN <- ifelse(covars$ctcae_grade >= 3, "Y", "N")
# # ------------------------------------------------------------------------------
# # 4. Could you please use the treatment/case definition parameters set for 5 and
# # 6 above but test for interaction (logit[y] ~ covariates + RV_burden +
# # RV_burden * treatment), changing the definitions of "covariates" based on how
# # you defined the treatment set in 5
# # ------------------------------------------------------------------------------

##############
## EUROPEAN ##
##############
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//diabetes/gene-based-analysis/rare_variant_analysis_EUR.Rdata")

## Run rare variant analysis; loop through genes

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

sum(chrALL$SKATO < 1e-3)
# 5

write.table(chrALL, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-EUR_chrALL_SKAT.pval", sep="\t", col.names=T, row.names=F, quote=F)

#############
## AFRICAN ##
#############
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//diabetes/gene-based-analysis/rare_variant_analysis_AFR.Rdata")

## Run rare variant analysis; loop through genes
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


write.table(chrALL, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-AFR_chrALL_SKAT.pval", sep="\t", col.names=T, row.names=F, quote=F)

