################################
## Create variables for Cindy ##
################################
# maxabdrtdose: Maximum dose to the abdomen + TBI, cGy
# maxpelvisrtdose: Maximum dose to the pelvis + TBI, cGy
# aa_class_dose_5: Cumulative Alkylating Agent: Classic (CED mg/m2) within 5 years of primary cancer diagnosis
# aa_class_dose_any: Cumulative Alkylating Agent: Classic (CED mg/m2)


# # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# # 5. Could you please re-evaluate the top findings from this RV analysis with genes with SKAT-O P<1e-3 in treatment-stratified samples, per the following:; Re: 5
# # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# # 5.a. Cutoff: for exposed, define as >200 cGy
# # 5.b. Cutoff: same as above, but for RT cutoffs at >=500 cGy vs <500, >=1000 cGy vs <1000, >=1500 cGy vs
# # 5.c Cutoff : same as above, but abdominal OR pelvic RT exposed

# # abdominalRT_YN 
# com.data$maxabdrtdose.exposed_more_than_200cGy_YN <- ifelse(com.data$maxabdrtdose > 200, 1, 0)
# com.data$maxabdrtdose.exposed_500cGy_or_higher_YN <- ifelse(com.data$maxabdrtdose >= 500, 1, 0)
# com.data$maxabdrtdose.exposed_1000cGy_or_higher_YN <- ifelse(com.data$maxabdrtdose >= 1000, 1, 0)
# com.data$maxabdrtdose.exposed_1500cGy_or_higher_YN <- ifelse(com.data$maxabdrtdose >= 1500, 1, 0)
# com.data$maxabdrtdose.exposed_2000cGy_or_higher_YN <- ifelse(com.data$maxabdrtdose >= 2000, 1, 0)
# 
# # PelvisRT_YN
# com.data$maxpelvisrtdose.exposed_more_than_200cGy_YN <- ifelse(com.data$maxpelvisrtdose > 200, 1, 0)
# com.data$maxpelvisrtdose.exposed_500cGy_or_higher_YN <- ifelse(com.data$maxpelvisrtdose >= 500, 1, 0)
# com.data$maxpelvisrtdose.exposed_1000cGy_or_higher_YN <- ifelse(com.data$maxpelvisrtdose >= 1000, 1, 0)
# com.data$maxpelvisrtdose.exposed_1500cGy_or_higher_YN <- ifelse(com.data$maxpelvisrtdose >= 1500, 1, 0)
# com.data$maxpelvisrtdose.exposed_2000cGy_or_higher_YN <- ifelse(com.data$maxpelvisrtdose >= 2000, 1, 0)
# 
# # Abdominal or Pelvic RT exposed
# com.data$abd_OR_pelvis_exposed_more_than_200cGY_YN <- ifelse(com.data$maxabdrtdose.exposed_more_than_200cGy_YN == 1 |com.data$maxpelvisrtdose.exposed_more_than_200cGy_YN == 1, 1,0)
# com.data$abd_OR_pelvis.exposed_500cGy_or_higher_YN <- ifelse(com.data$maxabdrtdose.exposed_500cGy_or_higher_YN == 1 |com.data$maxpelvisrtdose.exposed_500cGy_or_higher_YN == 1, 1,0)
# com.data$abd_OR_pelvis.exposed_1000cGy_or_higher_YN <- ifelse(com.data$maxabdrtdose.exposed_1000cGy_or_higher_YN == 1 |com.data$maxpelvisrtdose.exposed_1000cGy_or_higher_YN == 1, 1,0)
# com.data$abd_OR_pelvis.exposed_1500cGy_or_higher_YN <- ifelse(com.data$maxabdrtdose.exposed_1500cGy_or_higher_YN == 1 |com.data$maxpelvisrtdose.exposed_1500cGy_or_higher_YN == 1, 1,0)
# com.data$abd_OR_pelvis.exposed_2000cGy_or_higher_YN <- ifelse(com.data$maxabdrtdose.exposed_2000cGy_or_higher_YN == 1 |com.data$maxpelvisrtdose.exposed_2000cGy_or_higher_YN == 1, 1,0)
# 
# # ----------------------
# # 5.d. alkylating agents
# # ----------------------
# # AA within 5 years of primary cancer
# com.data$aa_class_dose_5 <- PHENO.ANY_SN$aa_class_dose_5[match(com.data$IID, PHENO.ANY_SN$sjlid)]
# com.data$aa_class_dose_5_YN <- ifelse(com.data$aa_class_dose_5 > 0, 1, 0 )
# com.data$aa_class_dose_5_4000_or_higher_YN <- ifelse(com.data$aa_class_dose_5 >= 4000, 1, 0 )
# # AA Any
# com.data$aa_class_dose_any <- PHENO.ANY_SN$aa_class_dose_any[match(com.data$IID, PHENO.ANY_SN$sjlid)]
# com.data$aa_class_dose_any_YN <- ifelse(com.data$aa_class_dose_any > 0, 1, 0 )
# com.data$aa_class_dose_any_4000_or_higher_YN <- ifelse(com.data$aa_class_dose_any >= 4000, 1, 0 )
# 
# # ------------------------------------------------------------------------------------
# # 6. Could you please repeat step 5 but with cases defined using CTCAE-graded DM  >=3?
# # ------------------------------------------------------------------------------------
# attach("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/CTCAE Abnormal Glucose Metabolism Rows With Transient Hyperglycemia-Mapping Rows Removed For Yadav's Population (n=4747, rows=9349).Rda")
# diabetes29.ctcae <- diabetes29.ctcae
# detach("file:Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/CTCAE Abnormal Glucose Metabolism Rows With Transient Hyperglycemia-Mapping Rows Removed For Yadav's Population (n=4747, rows=9349).Rda")
# dim(diabetes29.ctcae)
# 
# ## Extract samples with the most recent grade dates
# diabetes29.ctcae <- diabetes29.ctcae[tapply(1:nrow(diabetes29.ctcae),diabetes29.ctcae$sjlid,function(ii) ii[which.max(diabetes29.ctcae$grade.date[ii])]),]
# sum(com.data$IID %in% com.data$IID)
# # 3113
# 
# com.data$ctcae_grade <- diabetes29.ctcae$grade[match(com.data$IID, diabetes29.ctcae$sjlid)]
# com.data$ctcae_grade[com.data$ctcae_grade == -9 ] <- NA
# com.data$ctcae_grad_3_or_higher_YN <- ifelse(com.data$ctcae_grade >= 3, 1, 0)
# # ------------------------------------------------------------------------------
# # 7. Could you please use the treatment/case definition parameters set for 5 and
# # 6 above but test for interaction (logit[y] ~ covariates + RV_burden +
# # RV_burden * treatment), changing the definitions of "covariates" based on how
# # you defined the treatment set in 5
# # com.data$RV_Burden <- rowSums(SETgeno) # created for each gene below
# # ------------------------------------------------------------------------------

##############
## EUROPEAN ##
##############
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//diabetes/gene-based-analysis/rare_variant_analysis_EUR.Rdata")
#-----------------------------------------------------------------------
## Run rare variant analysis using covars as in GWAS; loop through genes
#-----------------------------------------------------------------------
## Initialize values:
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


chrALL.EUR <- cbind.data.frame(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
colnames(chrALL.EUR) <- chrALL.EUR[1,]
chrALL.EUR <- chrALL.EUR[-1,]
chrALL.EUR[,c("SKAT", "SKATO", "BURDEN")] <- sapply(chrALL.EUR[,c("SKAT", "SKATO", "BURDEN")], as.numeric)



# write.table(chrALL.EUR, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-EUR_chrALL_SKAT.pval", sep="\t", col.names=T, row.names=F, quote=F)


#################################################################
## Additional analysis on selected genes as requested by Cindy ##
#################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# 5 - Could you please re-evaluate the top findings from this RV analysis with genes with SKAT-O P<1e-3 in treatment-stratified samples, per the following:
#   abdominal RT exposed (for exposed, define as >200 cGy)
# pelvic RT exposed
# abdominal OR pelvic RT exposed
# same as above, but for RT cutoffs at >=500 cGy vs <500, >=1000 cGy vs <1000, >=1500 cGy vs <1500, >=2000 cGy vs <2000
# alkylating agents exposed (none vs any)
# alkylating agents, using dose cutoffs at >=4000 mg/m2 vs <4000 mg/m2
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# selecting genes with SKATO 1e-3
sum(chrALL.EUR$SKATO < 1e-3)
# 5


wanted.genes <- chrALL.EUR$GENE[chrALL.EUR$SKATO < 1e-3]
geneLIST.wanted <- geneLIST[names(geneLIST) %in% wanted.genes]

# Loop over different treatment cutoffs requested by Cindy for request #5 and #6
treatments <- c("maxabdrtdose.exposed_more_than_200cGy_YN", # abdominal RT
            "maxabdrtdose.exposed_500cGy_or_higher_YN",
            "maxabdrtdose.exposed_1000cGy_or_higher_YN",
            "maxabdrtdose.exposed_1500cGy_or_higher_YN",
            "maxabdrtdose.exposed_2000cGy_or_higher_YN",
            "aa_class_dose_5_YN", # aa_class_dose_5 (none vs any)
            "aa_class_dose_5_4000_or_higher_YN") # aa_class_dose_5 (>=4000 mg/m2 vs <4000 mg/m2)
 

### Exposed to treatments
chrALL.exposed <- NULL   # stratified by treatment (exposed)         
for (j in 1:length(treatments)){

# subset data
com.data.sub <- com.data[com.data[treatments[j]] == "1"  ,] # exposed to treatment
  
genelist<-"GENE"
SNPID<-"SNPID"
N0<-"N0"
N1<-"N1"
SKAT<-"SKAT"
SKATO<-"SKATO"
BURDEN<-"BURDEN"
### START loop over list of genes
for (i in 1:length(geneLIST.wanted)) {
  print(i)
  genename=names(geneLIST.wanted[i])
  print(genename)
  
  ### GET GENO DATA
  SETlist <- unname(rrapply(geneLIST.wanted[i],  how = 'unlist'))
  snpID<-paste(SETlist, collapse = ";")
  SETgeno<-as.matrix(com.data.sub[,SETlist])
  
  ### Define y
  y=as.matrix(com.data.sub$t2d)
  table(y)

  
  # Model adjust
  covariates = c("agedx","gender","age_last_visit","BMIadj","aa_class_dose_5","maxabdrtdose","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  

  index= which(colnames(com.data.sub)%in%covariates)
  Xm2=as.matrix(com.data.sub[,index])
  
  
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


chrALL.tmp <- cbind.data.frame(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
colnames(chrALL.tmp) <- chrALL.tmp[1,]
chrALL.tmp <- chrALL.tmp[-1,]
chrALL.tmp[,c("SKAT", "SKATO", "BURDEN")] <- sapply(chrALL.tmp[,c("SKAT", "SKATO", "BURDEN")], as.numeric)
chrALL.tmp$TREATMENT <- treatments[j]
chrALL.exposed <- rbind.data.frame(chrALL.exposed, chrALL.tmp)
}

write.table(chrALL.exposed, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-EUR_5.genes.exposed.to.treatments.t2d.as.caco.pval", sep="\t", col.names=T, row.names=F, quote=F)

### Not exposed to treatments
chrALL.not.exposed <- NULL   # stratified by treatment (not exposed)         
for (j in 1:length(treatments)){
  
  # subset data
  com.data.sub <- com.data[com.data[treatments[j]] == "0"  ,] # Not exposed to treatment
  
  genelist<-"GENE"
  SNPID<-"SNPID"
  N0<-"N0"
  N1<-"N1"
  SKAT<-"SKAT"
  SKATO<-"SKATO"
  BURDEN<-"BURDEN"
  ### START loop over list of genes
  for (i in 1:length(geneLIST.wanted)) {
    print(i)
    genename=names(geneLIST.wanted[i])
    print(genename)
    
    ### GET GENO DATA
    SETlist <- unname(rrapply(geneLIST.wanted[i],  how = 'unlist'))
    snpID<-paste(SETlist, collapse = ";")
    SETgeno<-as.matrix(com.data.sub[,SETlist])
    
    ### Define y
    y=as.matrix(com.data.sub$t2d)
    table(y)
    
    
    # Model adjust
    covariates = c("agedx","gender","age_last_visit","BMIadj","aa_class_dose_5","maxabdrtdose","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
    
    
    index= which(colnames(com.data.sub)%in%covariates)
    Xm2=as.matrix(com.data.sub[,index])
    
    
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
  
  
  chrALL.tmp <- cbind.data.frame(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
  colnames(chrALL.tmp) <- chrALL.tmp[1,]
  chrALL.tmp <- chrALL.tmp[-1,]
  chrALL.tmp[,c("SKAT", "SKATO", "BURDEN")] <- sapply(chrALL.tmp[,c("SKAT", "SKATO", "BURDEN")], as.numeric)
  chrALL.tmp$TREATMENT <- treatments[j]
  chrALL.not.exposed <- rbind.data.frame(chrALL.not.exposed, chrALL.tmp)
}

write.table(chrALL.not.exposed, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-EUR_5.genes.not.exposed.to.treatments.t2d.as.caco.pval", sep="\t", col.names=T, row.names=F, quote=F)

#--------------------------------------------------------------------------------------
# 6 - Could you please repeat step 5 but with cases defined using CTCAE-graded DM  >=3?
# -------------------------------------------------------------------------------------
### Exposed to treatments
chrALL.2.exposed <- NULL            
for (j in 1:length(treatments)){
  
  # subset data
  com.data.sub <- com.data[com.data[treatments[j]] == "0"  ,] # Not exposed to treatment
  
  
  genelist<-"GENE"
  SNPID<-"SNPID"
  N0<-"N0"
  N1<-"N1"
  SKAT<-"SKAT"
  SKATO<-"SKATO"
  BURDEN<-"BURDEN"
  ### START loop over list of genes
  for (i in 1:length(geneLIST.wanted)) {
    print(i)
    genename=names(geneLIST.wanted[i])
    print(genename)
    
    ### GET GENO DATA
    SETlist <- unname(rrapply(geneLIST.wanted[i],  how = 'unlist'))
    snpID<-paste(SETlist, collapse = ";")
    SETgeno<-as.matrix(com.data.sub[,SETlist])
    
    ### Define y
    # # 6 - Could you please repeat step 5 but with cases defined using CTCAE-graded DM  >=3?
    y=as.matrix(com.data.sub$ctcae_grad_3_or_higher_YN)
    table(y)
    
    # Model adjust
    covariates = c("agedx","gender","age_last_visit","BMIadj","aa_class_dose_5","maxabdrtdose","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
    
    
    index= which(colnames(com.data.sub)%in%covariates)
    Xm2=as.matrix(com.data.sub[,index])
    
    
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
  
  
  chrALL.tmp <- cbind.data.frame(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
  colnames(chrALL.tmp) <- chrALL.tmp[1,]
  chrALL.tmp <- chrALL.tmp[-1,]
  chrALL.tmp[,c("SKAT", "SKATO", "BURDEN")] <- sapply(chrALL.tmp[,c("SKAT", "SKATO", "BURDEN")], as.numeric)
  chrALL.tmp$TREATMENT <- treatments[j]
  chrALL.2.exposed <- rbind.data.frame(chrALL.2.exposed, chrALL.tmp)
}

write.table(chrALL.2.exposed, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-EUR_5.genes.exposed.to.treatments.CTCAE.gt.3.as.caco.pval", sep="\t", col.names=T, row.names=F, quote=F)


### Not exposed to treatments
chrALL.2.not.exposed <- NULL            
for (j in 1:length(treatments)){
  
  # subset data
  com.data.sub <- com.data[com.data[treatments[j]] == "1"  ,] # exposed to treatment
  
  
  genelist<-"GENE"
  SNPID<-"SNPID"
  N0<-"N0"
  N1<-"N1"
  SKAT<-"SKAT"
  SKATO<-"SKATO"
  BURDEN<-"BURDEN"
  ### START loop over list of genes
  for (i in 1:length(geneLIST.wanted)) {
    print(i)
    genename=names(geneLIST.wanted[i])
    print(genename)
    
    ### GET GENO DATA
    SETlist <- unname(rrapply(geneLIST.wanted[i],  how = 'unlist'))
    snpID<-paste(SETlist, collapse = ";")
    SETgeno<-as.matrix(com.data.sub[,SETlist])
    
    ### Define y
    # # 6 - Could you please repeat step 5 but with cases defined using CTCAE-graded DM  >=3?
    y=as.matrix(com.data.sub$ctcae_grad_3_or_higher_YN)
    table(y)
    
    # Model adjust
    covariates = c("agedx","gender","age_last_visit","BMIadj","aa_class_dose_5","maxabdrtdose","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
    
    
    index= which(colnames(com.data.sub)%in%covariates)
    Xm2=as.matrix(com.data.sub[,index])
    
    
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
  
  
  chrALL.tmp <- cbind.data.frame(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
  colnames(chrALL.tmp) <- chrALL.tmp[1,]
  chrALL.tmp <- chrALL.tmp[-1,]
  chrALL.tmp[,c("SKAT", "SKATO", "BURDEN")] <- sapply(chrALL.tmp[,c("SKAT", "SKATO", "BURDEN")], as.numeric)
  chrALL.tmp$TREATMENT <- treatments[j]
  chrALL.2.not.exposed <- rbind.data.frame(chrALL.2.not.exposed, chrALL.tmp)
}

write.table(chrALL.2.not.exposed, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-EUR_5.genes.not.exposed.to.treatments.CTCAE.gt.3.as.caco.pval", sep="\t", col.names=T, row.names=F, quote=F)

#-------------------------------------------------------------------------------
# 7. Could you please use the treatment/case definition parameters set for 5
# above but test for interaction (logit[y] ~ covariates + RV_burden +
# RV_burden * treatment), changing the definitions of "covariates" based on how
# you defined the treatment set in 5?
#-------------------------------------------------------------------------------
RV_BURDEN_INTERACTION <- "RV_BURDEN_INTERACTION"

### START loop over formula
formula = c("t2d ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_more_than_200cGy_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_more_than_200cGy_YN", # abdominal RT
            "t2d ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_500cGy_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_500cGy_or_higher_YN",
            "t2d ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_1000cGy_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_1000cGy_or_higher_YN",
            "t2d ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_1500cGy_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_1500cGy_or_higher_YN",
            "t2d ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_2000cGy_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_2000cGy_or_higher_YN",
            "t2d ~ agedx + gender + age_last_visit + BMIadj + aa_class_dose_5_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*aa_class_dose_5_YN", # aa_class_dose_5 (none vs any)
            "t2d ~ agedx + gender + age_last_visit + BMIadj + aa_class_dose_5_4000_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*aa_class_dose_5_4000_or_higher_YN") # aa_class_dose_5 (>=4000 mg/m2 vs <4000 mg/m2)
               
            
### Case/Control based on t2d variable

RV_BURDEN_INTERACTION <- NULL

for (h in 1:length (formula)){
### START loop over list of genes
for (i in 1:length(geneLIST.wanted)) {
  print(i)
  genename=names(geneLIST.wanted[i])
  print(genename)
  
  ### GET GENO DATA
  SETlist <- unname(rrapply(geneLIST.wanted[i],  how = 'unlist'))
  snpID<-paste(SETlist, collapse = ";")
  SETgeno<-as.matrix(com.data[,SETlist])
  
 
  com.data$RV_Burden <- rowSums(SETgeno)
  glm.fit <- glm(formula[h], family = binomial(link = "logit"), data = com.data)  
  
  fit.df <- setNames(as.data.frame(coef(summary(glm.fit))[,4]), "P")
  fit.df$VAR <- rownames(fit.df)
  
  fit.df <- fit.df[grepl("RV", fit.df$VAR),]
  
  for (j in 1:nrow(fit.df)){
  RV_BURDEN_INTERACTION.tmp <- cbind.data.frame(genename, fit.df[j,2],fit.df[j,1], formula[h])
  colnames(RV_BURDEN_INTERACTION.tmp) <- c("GENE", "VAR", "P", "formula")
  RV_BURDEN_INTERACTION <- rbind.data.frame(RV_BURDEN_INTERACTION, RV_BURDEN_INTERACTION.tmp)  
  }
  
}
}

write.table(RV_BURDEN_INTERACTION, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-EUR_5.genes.RV_BURDEN_INTERACTION.t2d.as.caco.pval", sep="\t", col.names=T, row.names=F, quote=F)

#-------------------------------------------------------------------------------
## 7. Could you please use the treatment/case definition parameters set for 
## 6 above but test for interaction (logit[y] ~ covariates + RV_burden +
## RV_burden * treatment), changing the definitions of "covariates" based on how
## you defined the treatment set in 5?
#-------------------------------------------------------------------------------

### Case/Control based on CTCAE-graded DM  >=3; START loop over formula
formula = c("ctcae_grad_3_or_higher_YN ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_more_than_200cGy_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_more_than_200cGy_YN", # abdominal RT
            "ctcae_grad_3_or_higher_YN ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_500cGy_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_500cGy_or_higher_YN",
            "ctcae_grad_3_or_higher_YN ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_1000cGy_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_1000cGy_or_higher_YN",
            "ctcae_grad_3_or_higher_YN ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_1500cGy_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_1500cGy_or_higher_YN",
            "ctcae_grad_3_or_higher_YN ~ agedx + gender + age_last_visit + BMIadj + maxabdrtdose.exposed_2000cGy_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*maxabdrtdose.exposed_2000cGy_or_higher_YN",
            "ctcae_grad_3_or_higher_YN ~ agedx + gender + age_last_visit + BMIadj + aa_class_dose_5_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*aa_class_dose_5_YN", # aa_class_dose_5 (none vs any)
            "ctcae_grad_3_or_higher_YN ~ agedx + gender + age_last_visit + BMIadj + aa_class_dose_5_4000_or_higher_YN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + RV_Burden + RV_Burden*aa_class_dose_5_4000_or_higher_YN") # aa_class_dose_5 (>=4000 mg/m2 vs <4000 mg/m2)


RV_BURDEN_INTERACTION.2 <- NULL

for (h in 1:length (formula)){
### START loop over list of genes
  for (i in 1:length(geneLIST.wanted)) {
    print(i)
    genename=names(geneLIST.wanted[i])
    print(genename)
    
    ### GET GENO DATA
    SETlist <- unname(rrapply(geneLIST.wanted[i],  how = 'unlist'))
    snpID<-paste(SETlist, collapse = ";")
    SETgeno<-as.matrix(com.data[,SETlist])
    
    
    com.data$RV_Burden <- rowSums(SETgeno)
    glm.fit <- glm(formula[h], family = binomial(link = "logit"), data = com.data)  
    
    fit.df <- setNames(as.data.frame(coef(summary(glm.fit))[,4]), "P")
    fit.df$VAR <- rownames(fit.df)
    
    fit.df <- fit.df[grepl("RV", fit.df$VAR),]
    
    for (j in 1:nrow(fit.df)){
      RV_BURDEN_INTERACTION.tmp <- cbind.data.frame(genename, fit.df[j,2],fit.df[j,1], formula[h])
      colnames(RV_BURDEN_INTERACTION.tmp) <- c("GENE", "VAR", "P", "formula")
      RV_BURDEN_INTERACTION.2 <- rbind.data.frame(RV_BURDEN_INTERACTION.2, RV_BURDEN_INTERACTION.tmp)  
    }
    
  }
}

write.table(RV_BURDEN_INTERACTION.2, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-EUR_5.genes.RV_BURDEN_INTERACTION.CTCAE.gt.3.as.caco.pval", sep="\t", col.names=T, row.names=F, quote=F)

#############
## AFRICAN ##
#############
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//diabetes/gene-based-analysis/rare_variant_analysis_AFR.Rdata")

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


chrALL.AFR <- cbind.data.frame(genelist,SNPID,N0,N1,SKAT,SKATO,BURDEN)
colnames(chrALL.AFR) <- chrALL.AFR[1,]
chrALL.AFR <- chrALL.AFR[-1,]
chrALL.AFR[,c("SKAT", "SKATO", "BURDEN")] <- sapply(chrALL.AFR[,c("SKAT", "SKATO", "BURDEN")], as.numeric)


# write.table(chrALL.AFR, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/gene-based-analysis/geneset-AFR_chrALL_SKAT.pval", sep="\t", col.names=T, row.names=F, quote=F)

