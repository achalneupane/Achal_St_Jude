# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.MENINGIOMA.V16.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "None"] <- "<30"
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- "<30"
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- "<30"
PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("<30", ">=30"))

table(PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$MENINGIOMA == 1])

PHENO.ANY_SN$epitxn_dose_5.category <- as.character(PHENO.ANY_SN$epitxn_dose_5.category)
PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "2nd"] <- "2nd-3rd"
PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "3rd"] <- "2nd-3rd"
PHENO.ANY_SN$epitxn_dose_5.category <- factor(PHENO.ANY_SN$epitxn_dose_5.category, levels = c("None", "1st", "2nd-3rd"))

#################
## Add new PRS ##
#################
testPRS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/Yadav_CH_related_06_28_2023/CH_related_PRSs_sjlife.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, testPRS[match(PHENO.ANY_SN$sjlid, testPRS$IID),])


newScores <- c("SCORE_chip", "SCORE_mCA", "SCORE_LOY", "SCORE_telomere_length", "SCORE_chip_cat", "SCORE_mCA_cat", "SCORE_LOY_cat", "SCORE_telomere_length_cat" )

PHENO.ANY_SN$SCORE_chip_cat <- factor(PHENO.ANY_SN$SCORE_chip_cat, levels = c("[-3.96,-0.662]", "(-0.662,0.701]", "(0.701,3]"))
PHENO.ANY_SN$SCORE_mCA_cat <- factor(PHENO.ANY_SN$SCORE_mCA_cat, levels = c("[-3.1,-0.69]", "(-0.69,0.687]", "(0.687,3.25]"))
PHENO.ANY_SN$SCORE_LOY_cat <- factor(PHENO.ANY_SN$SCORE_LOY_cat, levels = c("[-3.64,-0.676]", "(-0.676,0.678]", "(0.678,3.71]"))
PHENO.ANY_SN$SCORE_telomere_length_cat <- factor(PHENO.ANY_SN$SCORE_telomere_length_cat, levels = c("[-4.06,-0.67]", "(-0.67,0.644]", "(0.644,3.44]"))

######################################
## Attributable fraction of Any SNs ##
######################################
dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

final.model <- {}
for (i in 1:length(newScores)){
  print(paste0("Doing i: ", i))
  updatedFormula <- as.formula(paste("event ~", newScores[i], "+ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + 
                maxsegrtdose.category + epitxn_dose_5.category + 
                EAS + AFR +
                any_chemo_missing + any_rt_missing"))
  
  
  fit_all = glm(formula = updatedFormula,
                family = "poisson", offset = log(dat_all$PY), data = dat_all)
  
  summary(fit_all)
  
  
  (output <- summary(fit_all)$coefficients)
  as.data.frame(apply(output, 2, formatC, format="f", digits=4))
  # options(scipen=999)
  estimate <- format(round(output[,1],3), nsmall = 3)
  std.error <- format(round(output[,2],3), nsmall = 3)
  # P.val <- formatC(output[,4], format="G", digits=3)
  P.val <- output[,4]
  P.val[P.val < 0.001] <- "<0.001"
  P.val[!grepl("<", P.val)] <- format(round(as.numeric(P.val[!grepl("<", P.val)]), 3), nsmall = 3)
  sn.model <- (setNames(cbind.data.frame(estimate, std.error, P.val
  ), c("Estimate", "Std.error", "P")))
  sn.model <- sn.model[!grepl("AGE_AT_LAST_CONTACT", row.names(sn.model)),]
  sn.model
  empty_row <- data.frame(Variables = NA, Estimate = NA, Std.error = NA, P = NA)
  rownames(empty_row) <-  paste0("ROW_", newScores[i])
  # Add the empty row to the dataframe
  sn.model$Variables <- row.names(sn.model)
  row.names(sn.model) <- NULL
  sn.model <- rbind(empty_row, sn.model)
  # View(sn.model)
  
  final.model <- rbind(final.model, sn.model)
}