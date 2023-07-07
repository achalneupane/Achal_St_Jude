## Phenotype data for SJLIFE
rm(list=ls())
setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/')

## CTCAE grades from 20200430 data freeze
library(sas7bdat)
ctcae_20200430 = read.sas7bdat('/Volumes/data/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat')
# Select cardiomyopathy conditions
cmp = subset(ctcae_20200430, condition == "Cardiomyopathy")
# Those with grade 2 or higher
cmp.2plus = subset(cmp, grade >= 2)
# Get the maximum graded condition
cmp.2plus.sorted = cmp.2plus[order(cmp.2plus$sjlid, cmp.2plus$ageevent, decreasing = TRUE),]
cmp.2plus.uniq = cmp.2plus.sorted[!duplicated(cmp.2plus.sorted$sjlid),]
# Those with grade 0
cmp.0 = subset(cmp, grade == 0)
cmp.0.clean = subset(cmp.0, !(sjlid %in% cmp.2plus.uniq$sjlid))
cmp.0.clean.sorted = cmp.0.clean[order(cmp.0.clean$sjlid, cmp.0.clean$ageevent, decreasing = TRUE),] # the most recent age at last contact among multiple controls
cmp.0.clean.uniq = cmp.0.clean.sorted[!duplicated(cmp.0.clean.sorted$sjlid),]
cmp.final = rbind(cmp.2plus.uniq, cmp.0.clean.uniq)
cmp.final = cmp.final[c('sjlid', 'grade', 'ageevent', 'gradedt', 'ejection_fraction')]

## Get diagnosis and demographics
diag = read.sas7bdat('/Volumes/data/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat')
diag_primx = subset(diag, primdx==1, select=c('sjlid','diagdt','agedx','diaggrp'))
demo = read.sas7bdat('/Volumes/data/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat')
diag_primx_demo = merge(diag_primx, demo, by='sjlid')
agelcontact = read.sas7bdat('/Volumes/data/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Tracking Data/lstcondt.sas7bdat')
agelcontact = agelcontact[c('sjlid', 'agelstcontact')]
clinical = merge(diag_primx_demo, agelcontact, by='sjlid')

## Treatment data
chemo_dose = read.sas7bdat('/Volumes/data/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat')
chemo_dose_anthra = chemo_dose[c('sjlid','anthra_jco_dose_any')]
chemo_yn = read.sas7bdat('/Volumes/data/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_yn.sas7bdat')
chemo_yn_anthra = chemo_yn[c('sjlid','doxorubicin_any','daunorubicin_any')]
chemo_anthra = merge(chemo_dose_anthra, chemo_yn_anthra, by='sjlid')
rt = read.sas7bdat('/Volumes/data/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiationsum_yn.sas7bdat')
rt_chest = rt[c('sjlid', 'Chest')]
tx = merge(chemo_anthra, rt_chest, by='sjlid')

## Clinical and tx data
clinical_tx = merge(clinical, tx, by.x = 'sjlid', all.x = TRUE)

## Add cardiomyopathy data
clinical_tx_cmp = merge(clinical_tx, cmp.final, by='sjlid')
clinical_tx_cmp$CMP2plus = ifelse(clinical_tx_cmp$grade >=2, 2, 1)
clinical_tx_cmp$CMP3plus = ifelse(clinical_tx_cmp$grade >=3, 2, ifelse(clinical_tx_cmp$grade == 0, 1, NA))
clinical_tx_cmp$gender[clinical_tx_cmp$gender=="Male"] = 1
clinical_tx_cmp$gender[clinical_tx_cmp$gender=="Female"] = 2
  
## Ancestry proportion
ancestry =  read.table('/Volumes/clusterhome/ysapkota/Work/sjlife_1_n_2/ancestry/sjlife_wgs_admixture.txt', header = TRUE)
final = merge(clinical_tx_cmp, ancestry, by = 'sjlid')
# Save everything for future re-analyses
write.table(final, 'sjlife_all.txt', row.names=F, quote=F, sep="\t")

## Based on previous JNCI analyses, genetic effects were stronger among survivors exposed to doxorubicin only (no chest RT and daunorubicin)
final_dox = subset(final, doxorubicin_any==1 & Chest==0 & daunorubicin_any==0)

## European ancestry (given the dox only subset is smaller, use admixture eur prop of >=0.8 to define European ancestry - this likely provides more survivors than PCA based definition)
final_dox_eur = subset(final_dox, EUR >= 0.8)

## Subset including required variables only and format to PLINK to allow addition of top PCs only
plinkfile = final_dox_eur[c('sjlid', 'sjlid', 'CMP2plus', 'agedx','gender','agelstcontact', 'anthra_jco_dose_any', 'Chest')]
plinkfile = plinkfile[complete.cases(plinkfile),]
colnames(plinkfile)[c(1,2)] = c('FID', 'IID')
write.table(plinkfile, 'sjlife_eur_dox_only.txt', row.names = FALSE, quote = FALSE)

## African ancestry
final_dox_afr = subset(final_dox, AFR >= 0.8)
plinkfile_afr = final_dox_afr[c('sjlid', 'sjlid', 'CMP2plus', 'agedx','gender','agelstcontact', 'anthra_jco_dose_any', 'Chest')]
plinkfile_afr = plinkfile_afr[complete.cases(plinkfile_afr),]
colnames(plinkfile_afr)[c(1,2)] = c('FID', 'IID')
write.table(plinkfile_afr, 'sjlife_afr_dox_only.txt', row.names = FALSE, quote = FALSE)


