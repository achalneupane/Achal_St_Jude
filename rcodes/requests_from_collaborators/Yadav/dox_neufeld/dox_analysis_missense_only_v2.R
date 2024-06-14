# In V2: Can you please provide distribution of survivors who had 2 or more variants,
  # assess if that is associated with risk of cardiomyopathy? While performing the
  # analysis, you should compare survivors with 2 or more variants vs those with
  # no variants.

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/yadav_dox_abstract_06_03_2024/")
df <- read.table("dox_abstract.txt", header = T, sep = "\t")
df$ID <- sub(";.*", "", df$ID)

####################
## Read geno file ##
####################
library(data.table)
library(haven)
raw <- as.data.frame(fread("extract_vars_dox_recodeA.raw", header = T) )
rownames(raw) <- raw$IID
raw <- raw[,-c(1:6)]
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
colnames(raw) <- gsub("\\.", ":", HEADER)
sum(df$ID %in% colnames(raw))
# 67172

####################
## Non-Synonymous ##
####################
df.nonsynonymous = df[grepl('missense', df$ANN....EFFECT),]
# df.nonsynonymous = df[grepl('stop_gain|frame|missense', df$ANN....EFFECT),]
table(df.nonsynonymous$ANN....EFFECT)
table(df.nonsynonymous$ANN....GENE)
raw.nonsynonymous <- raw[colnames(raw) %in% df.nonsynonymous$ID]
raw.nonsynonymous$varSum <- rowSums(raw.nonsynonymous)
raw.nonsynonymous$carrier <- ifelse (raw.nonsynonymous$varSum > 0, 1,0)

df.nonsynonymous.unique <- df.nonsynonymous[!duplicated(df.nonsynonymous$ID),]
table(df.nonsynonymous.unique$ANN....GENE)
# CHPT1    GNPTAB   SLC25A3 UHRF1BP1L 
# 24        57        21        83 

##########
## P/LP ##
##########
df.clinvar <- df[grepl('^Pathogenic', df$CLNSIG),]
table(df.clinvar$ANN....EFFECT)
table(df.clinvar$ANN....GENE)

raw.clinvar <- raw[colnames(raw) %in% df.clinvar$ID]
raw.clinvar$varSum <- rowSums(raw.clinvar)
raw.clinvar$carrier <- ifelse (raw.clinvar$varSum > 0, 1,0)
#############################################
## Pheno with anthra > 0 and chestRT < 200 ##
#############################################
# pheno.sjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3.pheno", header = T)
pheno.sjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/Updated_CMP_phenotype_all_survivors.txt", header = T)
radiation <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiation_dosimetry.sas7bdat')
pheno.sjlife$chestRT <- radiation$maxchestrtdose[match(pheno.sjlife$IID, radiation$sjlid)]
pheno.sjlife <- pheno.sjlife[which(pheno.sjlife$chestRT <=200),]
pheno.sjlife <- pheno.sjlife[which(pheno.sjlife$anthra_jco_dose_any > 0),]

#####################################
## Add carrier status to phenotype ##
#####################################
## nonsynonymous
pheno.sjlife$carrier_nonsynonymous <- raw.nonsynonymous$carrier[match(pheno.sjlife$IID, rownames(raw.nonsynonymous))]
## clinvar
pheno.sjlife$carrier_PLP_clinvar <- raw.clinvar$carrier[match(pheno.sjlife$IID, rownames(raw.clinvar))]


#Can you please provide distribution of survivors who had 2 or more variants,
#assess if that is associated with risk of cardiomyopathy? While performing the
#analysis, you should compare survivors with 2 or more variants vs those with no
#variants. 

# nonsynonymous
pheno.sjlife$varSum.nonsynonymous <- raw.nonsynonymous$varSum[match(pheno.sjlife$IID, rownames(raw.nonsynonymous))]
table(pheno.sjlife$varSum.nonsynonymous)
## clinvar
pheno.sjlife$varSum.clinvar <- raw.clinvar$varSum[match(pheno.sjlife$IID, rownames(raw.clinvar))]


## With more than 2 carriers
pheno.sjlife <- pheno.sjlife[pheno.sjlife$varSum.nonsynonymous!=1,]
as.data.frame(table(pheno.sjlife$varSum.nonsynonymous))

pheno.sjlife$carrier_nonsynonymous <- ifelse(pheno.sjlife$varSum.nonsynonymous > 0, 1, 0)


###########################
## Fisher test CMP2 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_nonsynonymous)
# p-value = 0.8052
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7481815 1.4695316
# sample estimates:
#   odds ratio 
# 1.046361 
table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_nonsynonymous)
# 0    1
# 1  406 582
# 2   72  108

###########################
## Fisher test CMP3 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_nonsynonymous)
# p-value = 0.6927
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5226287 1.5800557
# sample estimates:
#   odds ratio 
# 0.9043776 
table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_nonsynonymous)
# 0    1
# 1  406 582
# 2   27   35

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# Yadav: Can you please add the number of carriers among survivors with and
# without cardiomyopathy as well? Also, could you describe the criteria to
# define these variants? In addition, can you perform the analyses stratified by
# ancestry?

##############
## EUR only ##
##############
EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/final_EUR-PCAS_EUR.eigenvec", header = T)
pheno.sjlife.EUR <- pheno.sjlife[pheno.sjlife$IID %in% EUR$IID,]

###########################
## Fisher test CMP2 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife.EUR$CMP2plus, pheno.sjlife.EUR$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife.EUR$CMP2plus, pheno.sjlife.EUR$carrier_nonsynonymous)
# p-value = 0.7837
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7308498 1.5516457
# sample estimates:
#   odds ratio 
# 1.06329 
table(pheno.sjlife.EUR$CMP2plus, pheno.sjlife.EUR$carrier_nonsynonymous)
# 0   1
# 1 345 412
# 2  63  80


###########################
## Fisher test CMP3 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_nonsynonymous)
# p-value = 0.4641
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4311158 1.4967753
# sample estimates:
#   odds ratio 
# 0.8041061 

table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_nonsynonymous)
# 0   1
# 1 345 412
# 2  25  24

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

##############
## AFR only ##
##############
AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/final_AFR-PCAS_AFR.eigenvec", header = T)
pheno.sjlife.AFR <- pheno.sjlife[pheno.sjlife$IID %in% AFR$IID,]

###########################
## Fisher test CMP2 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife.AFR$CMP2plus, pheno.sjlife.AFR$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife.AFR$CMP2plus, pheno.sjlife.AFR$carrier_nonsynonymous)
# p-value = 0.7682
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.3865001 8.2503269
# sample estimates:
#   odds ratio 
# 1.458818 
table(pheno.sjlife.AFR$CMP2plus, pheno.sjlife.AFR$carrier_nonsynonymous)
# 0   1
# 1  25 114
# 2   3  20

###########################
## Fisher test CMP3 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife.AFR$CMP3plus, pheno.sjlife.AFR$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife.AFR$CMP3plus, pheno.sjlife.AFR$carrier_nonsynonymous)
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.217274 80.893622
# sample estimates:
#   odds ratio 
# 1.74873 
table(pheno.sjlife.AFR$CMP3plus, pheno.sjlife.AFR$carrier_nonsynonymous)
#      0   1
# 1   25   114
# 2   1    8





