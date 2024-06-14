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
df.nonsynonymous = df[grepl('stop_gain|frame|missense', df$ANN....EFFECT),]
# df.nonsynonymous = df[grepl('stop_gain|frame|missense', df$ANN....EFFECT),]
table(df.nonsynonymous$ANN....EFFECT)
table(df.nonsynonymous$ANN....GENE)
raw.nonsynonymous <- raw[colnames(raw) %in% df.nonsynonymous$ID]
raw.nonsynonymous$varSum <- rowSums(raw.nonsynonymous)
raw.nonsynonymous$carrier <- ifelse (raw.nonsynonymous$varSum > 0, 1,0)


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

###########################
## Fisher test CMP2 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_nonsynonymous)
# p-value = 0.6341
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.6778188 1.2824965
# sample estimates:
# odds ratio 
#  0.9275264
table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_nonsynonymous)
#     0    1
# 1  334 1272
# 2   62  219

PLP.test <- fisher.test(table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_PLP_clinvar))
PLP.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_PLP_clinvar)
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.00000 6.24878
# sample estimates:
#   odds ratio 
# 0 

table(pheno.sjlife$CMP2plus, pheno.sjlife$carrier_PLP_clinvar)
#      0    1
# 1 1601    5
# 2  281    0

###########################
## Fisher test CMP3 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_nonsynonymous)
# p-value = 0.2972
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.467314 1.294236
# sample estimates:
#   odds ratio 
# 0.7659829 
table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_nonsynonymous)
#      0    1
# 1  334 1272
# 2   24   70

PLP.test <- fisher.test(table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_PLP_clinvar))
PLP.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_PLP_clinvar)
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.00000 18.81983
# sample estimates:
#   odds ratio 
# 0

table(pheno.sjlife$CMP3plus, pheno.sjlife$carrier_PLP_clinvar)
#      0    1
# 1 1601    5
# 2   94    0

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
# p-value = 0.6014
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6490135 1.3110223
# sample estimates:
#   odds ratio 
# 0.9170534 
table(pheno.sjlife.EUR$CMP2plus, pheno.sjlife.EUR$carrier_nonsynonymous)
#     0   1
# 1 281 954
# 2  53 165

PLP.test <- fisher.test(table(pheno.sjlife.EUR$CMP2plus, pheno.sjlife.EUR$carrier_PLP_clinvar))
PLP.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife.EUR$CMP2plus, pheno.sjlife.EUR$carrier_PLP_clinvar)
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.000000 6.197173
# sample estimates:
#   odds ratio 
# 0 

table(pheno.sjlife.EUR$CMP2plus, pheno.sjlife.EUR$carrier_PLP_clinvar)
#      0    1
# 1 1230    5
# 2  218    0

###########################
## Fisher test CMP3 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_nonsynonymous)
# p-value = 0.1087
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.3731323 1.1387733
# sample estimates:
#   odds ratio 
# 0.6428935 
table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_nonsynonymous)
#      0    1
# 1  281   954
# 2   22   48

PLP.test <- fisher.test(table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_PLP_clinvar))
PLP.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_PLP_clinvar)
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.00000 19.49179
# sample estimates:
#   odds ratio 
# 0 

table(pheno.sjlife.EUR$CMP3plus, pheno.sjlife.EUR$carrier_PLP_clinvar)
#      0    1
# 1 1230    5
# 2   70    0



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
# p-value = 0.7777
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4320515 8.4777175
# sample estimates:
#   odds ratio 
# 1.548016
table(pheno.sjlife.AFR$CMP2plus, pheno.sjlife.AFR$carrier_nonsynonymous)
#     0   1
# 1 23    178
# 2  3    36

# PLP.test <- fisher.test(table(pheno.sjlife.AFR$CMP2plus, pheno.sjlife.AFR$carrier_PLP_clinvar))
# PLP.test

expanded_table <- table(factor(pheno.sjlife.AFR$CMP2plus, levels = c(1,2)), factor(pheno.sjlife.AFR$carrier_PLP_clinvar, levels = c(0,1)))
# Perform Fisher's exact test
PLP.test <- fisher.test(expanded_table)
PLP.test


table(pheno.sjlife.AFR$CMP2plus, pheno.sjlife.AFR$carrier_PLP_clinvar)
#      0  
# 1    201
# 2    39

###########################
## Fisher test CMP3 plus ##
###########################
synonymous.test <- fisher.test(table(pheno.sjlife.AFR$CMP3plus, pheno.sjlife.AFR$carrier_nonsynonymous))
synonymous.test
# Fisher's Exact Test for Count Data
# data:  table(pheno.sjlife.AFR$CMP3plus, pheno.sjlife.AFR$carrier_nonsynonymous)
# p-value = 0.701
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.2921682 90.3885779
# sample estimates:
#   odds ratio 
# 2.062019
table(pheno.sjlife.AFR$CMP3plus, pheno.sjlife.AFR$carrier_nonsynonymous)
#      0   1
# 1   23   178
# 2   1    16

# PLP.test <- fisher.test(table(pheno.sjlife.AFR$CMP3plus, pheno.sjlife.AFR$carrier_PLP_clinvar))
# PLP.test

expanded_table <- table(factor(pheno.sjlife.AFR$CMP3plus, levels = c(1,2)), factor(pheno.sjlife.AFR$carrier_PLP_clinvar, levels = c(0,1)))
# Perform Fisher's exact test
PLP.test <- fisher.test(expanded_table)
PLP.test


table(pheno.sjlife.AFR$CMP3plus, pheno.sjlife.AFR$carrier_PLP_clinvar)
#     0    
# 1   201  
# 2   17   

