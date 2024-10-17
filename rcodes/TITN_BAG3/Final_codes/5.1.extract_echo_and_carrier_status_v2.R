rm(list=ls())
merged.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno", header = T)
sjlife.eur.dat <- merged.dat
sjlife.afr.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3_afr_kendrick.pheno", header = T)

# Keep sjlife and ccss_exp
merged.dat2 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
sjlife.eur.dat2 <- merged.dat2[merged.dat2$cohort==1| merged.dat2$cohort==3,]
sjlife.eur.dat <- sjlife.eur.dat[sjlife.eur.dat$IID %in% sjlife.eur.dat2$IID,]
dim(sjlife.eur.dat)

# sjlife.afr.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3_afr.pheno", header = T)


colnames(sjlife.eur.dat)
colnames(sjlife.afr.dat)
sjlife.eur.dat <- sjlife.eur.dat[c("FID", "IID", "CMP", "agedx", "agelstcontact", "gender", "anthra_jco_dose_any", "hrtavg",
                 "ejection_fraction_hrt", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "cohort_two")]
sjlife.eur.dat$cohort <- NA
sjlife.eur.dat$cohort[sjlife.eur.dat$cohort_two == 1] <- "sjlife_eur"
sjlife.eur.dat$cohort[sjlife.eur.dat$cohort_two == 3] <- "ccss_exp_eur"
sjlife.eur.dat$ancestry <- "EUR"
sjlife.afr.dat <- sjlife.afr.dat[c("FID", "IID", "CMP", "agedx", "agelstcontact", "gender", "anthra_jco_dose_any", "hrtavg",
                                   "ejection_fraction_hrt", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
sjlife.afr.dat$cohort <- "sjlife_afr"
sjlife.afr.dat$ancestry <- "AFR"

EUR.dat.PLP <- sjlife.eur.dat
AFR.dat.PLP <- sjlife.afr.dat


## Add PLP EUR
sjlife.BAG3.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/BAG3_EUR_vars_recodeA.raw", header = T)
sjlife.BAG3.PLP.EUR$carrier <- ifelse(rowSums(sjlife.BAG3.PLP.EUR[grepl("chr", colnames(sjlife.BAG3.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.DSP.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/DSP_EUR_vars_recodeA.raw", header = T)
sjlife.DSP.PLP.EUR$carrier <- ifelse(rowSums(sjlife.DSP.PLP.EUR[grepl("chr", colnames(sjlife.DSP.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.LMNA.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/LMNA_EUR_vars_recodeA.raw", header = T)
sjlife.LMNA.PLP.EUR$carrier <- ifelse(rowSums(sjlife.LMNA.PLP.EUR[grepl("chr", colnames(sjlife.LMNA.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.MYH7.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/MYH7_EUR_vars_recodeA.raw", header = T)
sjlife.MYH7.PLP.EUR$carrier <- ifelse(rowSums(sjlife.MYH7.PLP.EUR[grepl("chr", colnames(sjlife.MYH7.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.SCN5A.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/SCN5A_EUR_vars_recodeA.raw", header = T)
sjlife.SCN5A.PLP.EUR$carrier <- ifelse(rowSums(sjlife.SCN5A.PLP.EUR[grepl("chr", colnames(sjlife.SCN5A.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.TCAP.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TCAP_EUR_vars_recodeA.raw", header = T)
sjlife.TCAP.PLP.EUR$carrier <- ifelse(rowSums(sjlife.TCAP.PLP.EUR[grepl("chr", colnames(sjlife.TCAP.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.TNNC1.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TNNC1_EUR_vars_recodeA.raw", header = T)
sjlife.TNNC1.PLP.EUR$carrier <- ifelse(rowSums(sjlife.TNNC1.PLP.EUR[grepl("chr", colnames(sjlife.TNNC1.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.TNNT2.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TNNT2_EUR_vars_recodeA.raw", header = T)
sjlife.TNNT2.PLP.EUR$carrier <- ifelse(rowSums(sjlife.TNNT2.PLP.EUR[grepl("chr", colnames(sjlife.TNNT2.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.TTN_PSI.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TTN_PSI_EUR_vars_recodeA.raw", header = T)
sjlife.TTN_PSI.PLP.EUR$carrier <- ifelse(rowSums(sjlife.TTN_PSI.PLP.EUR[grepl("chr", colnames(sjlife.TTN_PSI.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.TTN_PSI_A_Band.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TTN_PSI_A_Band_EUR_vars_recodeA.raw", header = T)
sjlife.TTN_PSI_A_Band.PLP.EUR$carrier <- ifelse(rowSums(sjlife.TTN_PSI_A_Band.PLP.EUR[grepl("chr", colnames(sjlife.TTN_PSI_A_Band.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.BAG3.TTN_PSI.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/BAG3.TTN_PSI_EUR_vars_recodeA.raw", header = T)
sjlife.BAG3.TTN_PSI.PLP.EUR$carrier <- ifelse(rowSums(sjlife.BAG3.TTN_PSI.PLP.EUR[grepl("chr", colnames(sjlife.BAG3.TTN_PSI.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.BAG3.TTN_PSI_A_Band.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/BAG3.TTN_PSI_A_Band_EUR_vars_recodeA.raw", header = T)
sjlife.BAG3.TTN_PSI_A_Band.PLP.EUR$carrier <- ifelse(rowSums(sjlife.BAG3.TTN_PSI_A_Band.PLP.EUR[grepl("chr", colnames(sjlife.BAG3.TTN_PSI_A_Band.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.ALL_PLP_With_TTN_PSI.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/ALL_PLP_With_TTN_PSI_EUR_vars_recodeA.raw", header = T)
sjlife.ALL_PLP_With_TTN_PSI.PLP.EUR$carrier <- ifelse(rowSums(sjlife.ALL_PLP_With_TTN_PSI.PLP.EUR[grepl("chr", colnames(sjlife.ALL_PLP_With_TTN_PSI.PLP.EUR))], na.rm = T) > 0, 1, 0)

sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/ALL_PLP_With_TTN_PSI_A_Band_EUR_vars_recodeA.raw", header = T)
sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.EUR$carrier <- ifelse(rowSums(sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.EUR[grepl("chr", colnames(sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.EUR))], na.rm = T) > 0, 1, 0)



EUR.dat.PLP$BAG3.PLP.carrier <- sjlife.BAG3.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.BAG3.PLP.EUR$IID)]
EUR.dat.PLP$DSP.PLP.carrier <- sjlife.DSP.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.DSP.PLP.EUR$IID)]
EUR.dat.PLP$LMNA.PLP.carrier <- sjlife.LMNA.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.LMNA.PLP.EUR$IID)]
EUR.dat.PLP$MYH7.PLP.carrier <- sjlife.MYH7.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.MYH7.PLP.EUR$IID)]
EUR.dat.PLP$SCN5A.PLP.carrier <- sjlife.SCN5A.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.SCN5A.PLP.EUR$IID)]
EUR.dat.PLP$TCAP.PLP.carrier <- sjlife.TCAP.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.TCAP.PLP.EUR$IID)]
EUR.dat.PLP$TNNC1.PLP.carrier <- sjlife.TNNC1.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.TNNC1.PLP.EUR$IID)]
EUR.dat.PLP$TNNT2.PLP.carrier <- sjlife.TNNT2.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.TNNT2.PLP.EUR$IID)]
EUR.dat.PLP$TTN_PSI.PLP.carrier <- sjlife.TTN_PSI.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.TTN_PSI.PLP.EUR$IID)]
EUR.dat.PLP$TTN_PSI_A_Band.PLP.carrier <- sjlife.TTN_PSI_A_Band.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.TTN_PSI_A_Band.PLP.EUR$IID)]
EUR.dat.PLP$BAG3.TTN_PSI.PLP.carrier <- sjlife.BAG3.TTN_PSI.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.BAG3.TTN_PSI.PLP.EUR$IID)]
EUR.dat.PLP$BAG3.TTN_PSI_A_Band.PLP.carrier <- sjlife.BAG3.TTN_PSI_A_Band.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.BAG3.TTN_PSI_A_Band.PLP.EUR$IID)]
EUR.dat.PLP$ALL_PLP_With_TTN_PSI.PLP.carrier <- sjlife.ALL_PLP_With_TTN_PSI.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.ALL_PLP_With_TTN_PSI.PLP.EUR$IID)]
EUR.dat.PLP$ALL_PLP_With_TTN_PSI_A_Band.PLP.carrier <- sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.EUR$carrier[match(EUR.dat.PLP$IID, sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.EUR$IID)]


## Add P/LP AFR
## Add PLP AFR
sjlife.BAG3.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/BAG3_AFR_vars_recodeA.raw", header = T)
sjlife.BAG3.PLP.AFR$carrier <- ifelse(rowSums(sjlife.BAG3.PLP.AFR[grepl("chr", colnames(sjlife.BAG3.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.DSP.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/DSP_AFR_vars_recodeA.raw", header = T)
sjlife.DSP.PLP.AFR$carrier <- ifelse(rowSums(sjlife.DSP.PLP.AFR[grepl("chr", colnames(sjlife.DSP.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.LMNA.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/LMNA_AFR_vars_recodeA.raw", header = T)
sjlife.LMNA.PLP.AFR$carrier <- ifelse(rowSums(sjlife.LMNA.PLP.AFR[grepl("chr", colnames(sjlife.LMNA.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.MYH7.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/MYH7_AFR_vars_recodeA.raw", header = T)
sjlife.MYH7.PLP.AFR$carrier <- ifelse(rowSums(sjlife.MYH7.PLP.AFR[grepl("chr", colnames(sjlife.MYH7.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.SCN5A.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/SCN5A_AFR_vars_recodeA.raw", header = T)
sjlife.SCN5A.PLP.AFR$carrier <- ifelse(rowSums(sjlife.SCN5A.PLP.AFR[grepl("chr", colnames(sjlife.SCN5A.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.TCAP.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TCAP_AFR_vars_recodeA.raw", header = T)
sjlife.TCAP.PLP.AFR$carrier <- ifelse(rowSums(sjlife.TCAP.PLP.AFR[grepl("chr", colnames(sjlife.TCAP.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.TNNC1.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TNNC1_AFR_vars_recodeA.raw", header = T)
sjlife.TNNC1.PLP.AFR$carrier <- ifelse(rowSums(sjlife.TNNC1.PLP.AFR[grepl("chr", colnames(sjlife.TNNC1.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.TNNT2.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TNNT2_AFR_vars_recodeA.raw", header = T)
sjlife.TNNT2.PLP.AFR$carrier <- ifelse(rowSums(sjlife.TNNT2.PLP.AFR[grepl("chr", colnames(sjlife.TNNT2.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.TTN_PSI.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TTN_PSI_AFR_vars_recodeA.raw", header = T)
sjlife.TTN_PSI.PLP.AFR$carrier <- ifelse(rowSums(sjlife.TTN_PSI.PLP.AFR[grepl("chr", colnames(sjlife.TTN_PSI.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.TTN_PSI_A_Band.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/TTN_PSI_A_Band_AFR_vars_recodeA.raw", header = T)
sjlife.TTN_PSI_A_Band.PLP.AFR$carrier <- ifelse(rowSums(sjlife.TTN_PSI_A_Band.PLP.AFR[grepl("chr", colnames(sjlife.TTN_PSI_A_Band.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.BAG3.TTN_PSI.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/BAG3.TTN_PSI_AFR_vars_recodeA.raw", header = T)
sjlife.BAG3.TTN_PSI.PLP.AFR$carrier <- ifelse(rowSums(sjlife.BAG3.TTN_PSI.PLP.AFR[grepl("chr", colnames(sjlife.BAG3.TTN_PSI.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.BAG3.TTN_PSI_A_Band.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/BAG3.TTN_PSI_A_Band_AFR_vars_recodeA.raw", header = T)
sjlife.BAG3.TTN_PSI_A_Band.PLP.AFR$carrier <- ifelse(rowSums(sjlife.BAG3.TTN_PSI_A_Band.PLP.AFR[grepl("chr", colnames(sjlife.BAG3.TTN_PSI_A_Band.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.ALL_PLP_With_TTN_PSI.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/ALL_PLP_With_TTN_PSI_AFR_vars_recodeA.raw", header = T)
sjlife.ALL_PLP_With_TTN_PSI.PLP.AFR$carrier <- ifelse(rowSums(sjlife.ALL_PLP_With_TTN_PSI.PLP.AFR[grepl("chr", colnames(sjlife.ALL_PLP_With_TTN_PSI.PLP.AFR))], na.rm = T) > 0, 1, 0)

sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/ALL_PLP_With_TTN_PSI_A_Band_AFR_vars_recodeA.raw", header = T)
sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.AFR$carrier <- ifelse(rowSums(sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.AFR[grepl("chr", colnames(sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.AFR))], na.rm = T) > 0, 1, 0)


AFR.dat.PLP$BAG3.PLP.carrier <- sjlife.BAG3.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.BAG3.PLP.AFR$IID)]
AFR.dat.PLP$DSP.PLP.carrier <- sjlife.DSP.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.DSP.PLP.AFR$IID)]
AFR.dat.PLP$LMNA.PLP.carrier <- sjlife.LMNA.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.LMNA.PLP.AFR$IID)]
AFR.dat.PLP$MYH7.PLP.carrier <- sjlife.MYH7.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.MYH7.PLP.AFR$IID)]
AFR.dat.PLP$SCN5A.PLP.carrier <- sjlife.SCN5A.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.SCN5A.PLP.AFR$IID)]
AFR.dat.PLP$TCAP.PLP.carrier <- sjlife.TCAP.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.TCAP.PLP.AFR$IID)]
AFR.dat.PLP$TNNC1.PLP.carrier <- sjlife.TNNC1.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.TNNC1.PLP.AFR$IID)]
AFR.dat.PLP$TNNT2.PLP.carrier <- sjlife.TNNT2.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.TNNT2.PLP.AFR$IID)]
AFR.dat.PLP$TTN_PSI.PLP.carrier <- sjlife.TTN_PSI.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.TTN_PSI.PLP.AFR$IID)]
AFR.dat.PLP$TTN_PSI_A_Band.PLP.carrier <- sjlife.TTN_PSI_A_Band.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.TTN_PSI_A_Band.PLP.AFR$IID)]
AFR.dat.PLP$BAG3.TTN_PSI.PLP.carrier <- sjlife.BAG3.TTN_PSI.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.BAG3.TTN_PSI.PLP.AFR$IID)]
AFR.dat.PLP$BAG3.TTN_PSI_A_Band.PLP.carrier <- sjlife.BAG3.TTN_PSI_A_Band.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.BAG3.TTN_PSI_A_Band.PLP.AFR$IID)]
AFR.dat.PLP$ALL_PLP_With_TTN_PSI.PLP.carrier <- sjlife.ALL_PLP_With_TTN_PSI.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.ALL_PLP_With_TTN_PSI.PLP.AFR$IID)]
AFR.dat.PLP$ALL_PLP_With_TTN_PSI_A_Band.PLP.carrier <- sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.AFR$carrier[match(AFR.dat.PLP$IID, sjlife.ALL_PLP_With_TTN_PSI_A_Band.PLP.AFR$IID)]


save(EUR.dat.PLP, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/EUR.dat.PLP.RData")
save(AFR.dat.PLP, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/AFR.dat.PLP.RData")

#######################################
## Extract common variant carriers
pheno_gwas = read.table('Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno', header = TRUE)
# Using sjlife_ccss_org_ccss_exp_samples_updated.bim

pheno_gwas$cohort[pheno_gwas$cohort == 1] <- "sjlife"
pheno_gwas$cohort[pheno_gwas$cohort == 2] <- "ccss_org"
pheno_gwas$cohort[pheno_gwas$cohort == 3] <- "ccss_exp"

raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/common_missense_vars_recodeA.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
# HEADER <- colnames(raw)[-c(1:6)]
HEADER=colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER

map.rsid <- read.table(text="Positions rsid
178562809 rs3829746
178633315	rs6723526
178580212	rs2303838
178566270	rs3731746
178571293	rs744426
178586693	rs2042996
178556967	rs9808377
178599800	rs1001238
178532834	rs3829747
178541464	rs3731749
178592420	rs16866406
178693639	rs2042995
178593864	rs2288569
178718769	rs16866465
178722403	rs12693166
178714366	rs12693164
178717600	rs13390491
178785681	rs35813871
178795185	rs16866538
178717810	rs2627043
178567458	rs12463674
178681132	rs36051007
178780128	rs10497520
178689578	rs2244492
178759031	rs2291310
178764734	rs2291311
178756224	rs7585334
178751160	rs922984
178741811	rs2627037
178663651	rs2163008
178747656	rs72648907
178710784	rs72648998
119676774	rs3858340 
119670121 rs2234962", header = T)

map.rsid$SNP <- NA
for (i in 1:nrow(map.rsid)){
map.rsid$SNP[i]   <-   colnames(raw)[grepl(map.rsid$Positions[i], colnames(raw))]
}

colnames(raw) <- map.rsid$rsid[match(colnames(raw), map.rsid$SNP)]
raw <- cbind.data.frame(ID=rownames(raw), raw)
rownames(raw) <- NULL

EUR_common_variants <- cbind.data.frame(pheno_gwas, raw[match(pheno_gwas$IID, raw$ID),])


## AFR
pheno_gwas <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3_afr_kendrick.pheno", header = T)
raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/common_missense_vars_AFR_recodeA.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
# HEADER <- colnames(raw)[-c(1:6)]
HEADER=colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER



colnames(raw) <- map.rsid$rsid[match(colnames(raw), map.rsid$SNP)]
raw <- cbind.data.frame(ID=rownames(raw), raw)
rownames(raw) <- NULL

AFR_common_variants <- cbind.data.frame(pheno_gwas, raw[match(pheno_gwas$IID, raw$ID),])


EUR_common_variants <- EUR_common_variants[!colnames(EUR_common_variants) %in% "ID"]
AFR_common_variants <- AFR_common_variants[!colnames(AFR_common_variants) %in% "ID"]

save(EUR_common_variants, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/EUR_common_variants.RData")
save(AFR_common_variants, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/AFR_common_variants.RData")
