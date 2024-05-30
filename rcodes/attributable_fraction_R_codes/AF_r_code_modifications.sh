## only for 18 d

find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$max/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$epitxn_dose_5.category/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$anthra_jco_dose_5.category/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$aa_class_dose_5.category/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$any_tx_missing/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$any_chemo_missing/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$any_rt_missing/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$any_lifestyle_missing/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$Current_smoker_yn/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$PhysicalActivity_yn/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$RiskyHeavyDrink_yn/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/^PHENO.ANY_SN\$Obese_yn/s/^/# /' {} \;

find . -type f -name '4*' -exec sed -i '/"any_tx_missing",/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/"any_rt_missing",/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/"any_chemo_missing",/s/^/# /' {} \;
find . -type f -name '4*' -exec sed -i '/"any_lifestyle_missing",/s/^/# /' {} \;

find . -type f -name '4*' -exec sed -i '/^# # #.*,$/d' {} \;

find . -type f -name '*model_fit_*' -exec sed -i 's/EAS + AFR +/EAS + AFR,/' {} \;


find . -type f -exec sed -i 's/18b/18d/' {} \;


add a new line after line dat_all = PHENO.ANY_SN. The line I want to insert is PHENO.ANY_SN[PHENO.ANY_SN$gender == "Female",]
find . -type f -name '*model_fit_BREASTcancer*' -exec sed -i '/dat_all = PHENO.ANY_SN/a dat_all = dat_all[dat_all$gender == "Female",]' {} \;



comment out line with $any_chemo_missing or $any_rt_missing or $any_lifestyle_missing
find . -type f -name '*model_fit*' -exec sed -i '/\$any_chemo_missing\|\$any_rt_missing\|\$any_lifestyle_missing/s/^/# /' {} \;


# rt
find . -type f -name '*model_fit*' -exec sed -i 's/dat_rt$maxsegrtdose\.category =/dat_rt$maxsegrtdose.category [!grepl("Unknown", dat_rt$maxsegrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_rt$maxabdrtdose\.category =/dat_rt$maxabdrtdose.category [!grepl("Unknown", dat_rt$maxabdrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_rt$maxchestrtdose\.category =/dat_rt$maxchestrtdose.category [!grepl("Unknown", dat_rt$maxchestrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_rt$maxneckrtdose\.category =/dat_rt$maxneckrtdose.category [!grepl("Unknown", dat_rt$maxneckrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_rt$maxpelvisrtdose\.category =/dat_rt$maxpelvisrtdose.category [!grepl("Unknown", dat_rt$maxpelvisrtdose.category)] =/g' {} \;

# tx
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx\$epitxn_dose_5\.category =/dat_tx\$epitxn_dose_5.category [!grepl("Unknown", dat_tx\$epitxn_dose_5.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx\$anthra_jco_dose_5\.category =/dat_tx\$anthra_jco_dose_5.category [!grepl("Unknown", dat_tx\$anthra_jco_dose_5.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx\$aa_class_dose_5\.category =/dat_tx\$aa_class_dose_5.category [!grepl("Unknown", dat_tx\$aa_class_dose_5.category)] =/g' {} \;


# trt and rt
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.rt$maxsegrtdose\.category =/dat_tx.rt$maxsegrtdose.category [!grepl("Unknown", dat_tx.rt$maxsegrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.rt$maxabdrtdose\.category =/dat_tx.rt$maxabdrtdose.category [!grepl("Unknown", dat_tx.rt$maxabdrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.rt$maxchestrtdose\.category =/dat_tx.rt$maxchestrtdose.category [!grepl("Unknown", dat_tx.rt$maxchestrtdose.category)] =/g' {} \;

find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.rt$maxneckrtdose\.category =/dat_tx.rt$maxneckrtdose.category [!grepl("Unknown", dat_tx.rt$maxneckrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.rt$maxpelvisrtdose\.category =/dat_tx.rt$maxpelvisrtdose.category [!grepl("Unknown", dat_tx.rt$maxpelvisrtdose.category)] =/g' {} \;

find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.rt\$epitxn_dose_5\.category =/dat_tx.rt\$epitxn_dose_5.category [!grepl("Unknown", dat_tx.rt\$epitxn_dose_5.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.rt\$anthra_jco_dose_5\.category =/dat_tx.rt\$anthra_jco_dose_5.category [!grepl("Unknown", dat_tx.rt\$anthra_jco_dose_5.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.rt\$aa_class_dose_5\.category =/dat_tx.rt\$aa_class_dose_5.category [!grepl("Unknown", dat_tx.rt\$aa_class_dose_5.category)] =/g' {} \;



## Lifestyle
find . -type f -name '*model_fit*' -exec sed -i 's/dat_lifestyle\$Current_smoker_yn = "No"/dat_lifestyle\$Current_smoker_yn [!grepl("Unknown", dat_lifestyle\$Current_smoker_yn)] = "No"/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_lifestyle\$PhysicalActivity_yn = "Yes"/dat_lifestyle\$PhysicalActivity_yn [!grepl("Unknown", dat_lifestyle\$PhysicalActivity_yn)] = "Yes"/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_lifestyle\$RiskyHeavyDrink_yn = "No"/dat_lifestyle\$RiskyHeavyDrink_yn [!grepl("Unknown", dat_lifestyle\$RiskyHeavyDrink_yn)] = "No"/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_lifestyle\$HEALTHY_Diet_yn = "Yes"/dat_lifestyle\$HEALTHY_Diet_yn [!grepl("Unknown", dat_lifestyle\$HEALTHY_Diet_yn)] = "Yes"/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_lifestyle\$Obese_yn = "No"/dat_lifestyle\$Obese_yn [!grepl("Unknown", dat_lifestyle\$Obese_yn)] = "No"/g' {} \;

###############################################################
## Now in combined
# rt
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxsegrtdose\.category =/dat_tx.prs.lifestyle$maxsegrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxsegrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxabdrtdose\.category =/dat_tx.prs.lifestyle$maxabdrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxabdrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxchestrtdose\.category =/dat_tx.prs.lifestyle$maxchestrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxchestrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxneckrtdose\.category =/dat_tx.prs.lifestyle$maxneckrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxneckrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxpelvisrtdose\.category =/dat_tx.prs.lifestyle$maxpelvisrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxpelvisrtdose.category)] =/g' {} \;

# tx
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx\$epitxn_dose_5\.category =/dat_tx\$epitxn_dose_5.category [!grepl("Unknown", dat_tx\$epitxn_dose_5.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx\$anthra_jco_dose_5\.category =/dat_tx\$anthra_jco_dose_5.category [!grepl("Unknown", dat_tx\$anthra_jco_dose_5.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx\$aa_class_dose_5\.category =/dat_tx\$aa_class_dose_5.category [!grepl("Unknown", dat_tx\$aa_class_dose_5.category)] =/g' {} \;


# trt and rt
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxsegrtdose\.category =/dat_tx.prs.lifestyle$maxsegrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxsegrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxabdrtdose\.category =/dat_tx.prs.lifestyle$maxabdrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxabdrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxchestrtdose\.category =/dat_tx.prs.lifestyle$maxchestrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxchestrtdose.category)] =/g' {} \;

find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxneckrtdose\.category =/dat_tx.prs.lifestyle$maxneckrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxneckrtdose.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle$maxpelvisrtdose\.category =/dat_tx.prs.lifestyle$maxpelvisrtdose.category [!grepl("Unknown", dat_tx.prs.lifestyle$maxpelvisrtdose.category)] =/g' {} \;

find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle\$epitxn_dose_5\.category =/dat_tx.prs.lifestyle\$epitxn_dose_5.category [!grepl("Unknown", dat_tx.prs.lifestyle\$epitxn_dose_5.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle\$anthra_jco_dose_5\.category =/dat_tx.prs.lifestyle\$anthra_jco_dose_5.category [!grepl("Unknown", dat_tx.prs.lifestyle\$anthra_jco_dose_5.category)] =/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle\$aa_class_dose_5\.category =/dat_tx.prs.lifestyle\$aa_class_dose_5.category [!grepl("Unknown", dat_tx.prs.lifestyle\$aa_class_dose_5.category)] =/g' {} \;



## Lifestyle
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle\$Current_smoker_yn = "No"/dat_tx.prs.lifestyle\$Current_smoker_yn [!grepl("Unknown", dat_tx.prs.lifestyle\$Current_smoker_yn)] = "No"/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle\$PhysicalActivity_yn = "Yes"/dat_tx.prs.lifestyle\$PhysicalActivity_yn [!grepl("Unknown", dat_tx.prs.lifestyle\$PhysicalActivity_yn)] = "Yes"/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle\$RiskyHeavyDrink_yn = "No"/dat_tx.prs.lifestyle\$RiskyHeavyDrink_yn [!grepl("Unknown", dat_tx.prs.lifestyle\$RiskyHeavyDrink_yn)] = "No"/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle\$HEALTHY_Diet_yn = "Yes"/dat_tx.prs.lifestyle\$HEALTHY_Diet_yn [!grepl("Unknown", dat_tx.prs.lifestyle\$HEALTHY_Diet_yn)] = "Yes"/g' {} \;
find . -type f -name '*model_fit*' -exec sed -i 's/dat_tx.prs.lifestyle\$Obese_yn = "No"/dat_tx.prs.lifestyle\$Obese_yn [!grepl("Unknown", dat_tx.prs.lifestyle\$Obese_yn)] = "No"/g' {} \;

#########################################################
## Replace 18b to add startified analysis for ancestry ##
#########################################################
# the final files will be saved as V19b
In these files:
find . -type f -name '*model_fit*', 
after this:
'^af_by_no_favorable_lifestyle.category.gteq.35$'

add the following lines:

'## EUR
N_no_favorable_lifestyle.category = sum(dat_all$admixture == "EUR"], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.EUR = (N_all.EUR - N_no_favorable_lifestyle.category) / N_all.EUR
af_by_no_favorable_lifestyle.category.EUR <- round(af_by_no_favorable_lifestyle.category.EUR,3)
af_by_no_favorable_lifestyle.category.EUR

## AFR
N_no_favorable_lifestyle.category = sum(dat_all$admixture == "AFR"], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.AFR = (N_all.AFR - N_no_favorable_lifestyle.category) / N_all.AFR
af_by_no_favorable_lifestyle.category.AFR <- round(af_by_no_favorable_lifestyle.category.AFR,3)
af_by_no_favorable_lifestyle.category.AFR'






## Add EUR and AFR variables
find . -type f -name '*model_fit*' -exec sed -i '/^filtered_cc[[:space:]]*$/a \\n## Admixture classification\nadmixture <- read.table("Z:\/ResearchHome\/Groups\/sapkogrp\/projects\/Genomics\/common\/\/sjlife\/MERGED_SJLIFE_1_2\/MERGED_SJLIFE_PLINK_PER_CHR\/PCA\/SJLIFE_4481_Admixture_PCA_ethnicity.csv", sep = "\\t", header = T)\nEUR.admix <- admixture$INDIVIDUAL[admixture$EUR > 0.8]\nAFR.admix <- admixture$INDIVIDUAL[admixture$AFR > 0.6]\n\nPHENO.ANY_SN$admixture <- NA\nPHENO.ANY_SN$admixture [PHENO.ANY_SN$sjlid %in% EUR.admix] <- "EUR"\nPHENO.ANY_SN$admixture [PHENO.ANY_SN$sjlid %in% AFR.admix] <- "AFR"' {} +

## add subset groups
find . -type f -name '*model_fit*' -exec sed -i '/N_all.gteq.35 = sum(dat_all$pred_all\[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35\], na.rm = TRUE) # subset by age 35/ a\## Subset by ancestry\nN_all.EUR = sum(dat_all$pred_all\[dat_all$admixture == "EUR"\], na.rm = TRUE) # subset by ancestry\nN_all.AFR = sum(dat_all$pred_all\[dat_all$admixture == "AFR"\], na.rm = TRUE) # subset by ancestry' {} +

## ADD for Tx
find . -type f -name '*model_fit*' -exec sed -i '/^[[:space:]]*af_by_tx\.gteq\.35[[:space:]]*$/a \\n## EUR\nN_no_tx = sum(dat_all$pred_no_tx[dat_all$admixture == "EUR"], na.rm = TRUE)\naf_by_tx.EUR = (N_all.EUR - N_no_tx) / N_all.EUR\naf_by_tx.EUR <- round(af_by_tx.EUR,3)\naf_by_tx.EUR\n\n## AFR\nN_no_tx = sum(dat_all$pred_no_tx[dat_all$admixture == "AFR"], na.rm = TRUE)\naf_by_tx.AFR = (N_all.AFR - N_no_tx) / N_all.AFR\naf_by_tx.AFR <- round(af_by_tx.AFR,3)\naf_by_tx.AFR' {} +


## Add for rt
find . -type f -name '*model_fit*' -exec sed -i '/^[[:space:]]*af_by_rt\.gteq\.35[[:space:]]*$/a \\n## EUR\nN_no_rt = sum(dat_all$pred_no_rt[dat_all$admixture == "EUR"], na.rm = TRUE)\naf_by_rt.EUR = (N_all.EUR - N_no_rt) / N_all.EUR\naf_by_rt.EUR <- round(af_by_rt.EUR,3)\naf_by_rt.EUR\n\n## AFR\nN_no_rt = sum(dat_all$pred_no_rt[dat_all$admixture == "AFR"], na.rm = TRUE)\naf_by_rt.AFR = (N_all.AFR - N_no_rt) / N_all.AFR\naf_by_rt.AFR <- round(af_by_rt.AFR,3)\naf_by_rt.AFR' {} +


## Add for missing RT
find . -type f -name '*model_fit*' -exec sed -i '/^af_by_rt\.gteq\.35 <- "-"/a \\n## EUR\naf_by_rt.EUR <- "-"\n\n## AFR\naf_by_rt.AFR <- "-"' {} +

## Add for tx and rt
find . -type f -name '*model_fit*' -exec sed -i '/^[[:space:]]*af_by_tx\.rt\.gteq\.35[[:space:]]*$/a \\n## EUR\nN_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$admixture == "EUR"], na.rm = TRUE)\naf_by_tx.rt.EUR = (N_all.EUR - N_no_tx.rt) / N_all.EUR\naf_by_tx.rt.EUR <- round(af_by_tx.rt.EUR,3)\naf_by_tx.rt.EUR\n\n## AFR\nN_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$admixture == "AFR"], na.rm = TRUE)\naf_by_tx.rt.AFR = (N_all.AFR - N_no_tx.rt) / N_all.AFR\naf_by_tx.rt.AFR <- round(af_by_tx.rt.AFR,3)\naf_by_tx.rt.AFR' {} +

## Add for PRS
find . -type f -name '*model_fit*' -exec sed -i '/^[[:space:]]*af_by_prs\.gteq\.35[[:space:]]*$/a \\n## EUR\nN_no_prs = sum(dat_all$pred_no_prs[dat_all$admixture == "EUR"], na.rm = TRUE)\naf_by_prs.EUR = (N_all.EUR - N_no_prs) / N_all.EUR\naf_by_prs.EUR <- round(af_by_prs.EUR,3)\naf_by_prs.EUR\n\n## AFR\nN_no_prs = sum(dat_all$pred_no_prs[dat_all$admixture == "AFR"], na.rm = TRUE)\naf_by_prs.AFR = (N_all.AFR - N_no_prs) / N_all.AFR\naf_by_prs.AFR <- round(af_by_prs.AFR,3)\naf_by_prs.AFR' {} +


## Add for combined
find . -type f -name '*model_fit*' -exec sed -i '/^[[:space:]]*af_by_combined\.gteq\.35[[:space:]]*$/a \\n## EUR\nN_no_combined = sum(dat_all$pred_no_combined[dat_all$admixture == "EUR"], na.rm = TRUE)\naf_by_combined.EUR = (N_all.EUR - N_no_combined) / N_all.EUR\naf_by_combined.EUR <- round(af_by_combined.EUR,3)\naf_by_combined.EUR\n\n## AFR\nN_no_combined = sum(dat_all$pred_no_combined[dat_all$admixture == "AFR"], na.rm = TRUE)\naf_by_combined.AFR = (N_all.AFR - N_no_combined) / N_all.AFR\naf_by_combined.AFR <- round(af_by_combined.AFR,3)\naf_by_combined.AFR' {} +


## Add lifestyle if no lifestle
find . -type f -name '*model_fit*' -exec sed -i '/^af_by_no_favorable_lifestyle\.category\.gteq\.35 <- "-"/a \\n## EUR\naf_by_no_favorable_lifestyle.category.EUR <- "-"\n\n## AFR\naf_by_no_favorable_lifestyle.category.AFR <- "-"' {} +


## ADD lifestyle if lifestles
find . -type f -name '*model_fit*' -exec sed -i '/^[[:space:]]*af_by_no_favorable_lifestyle\.category\.gteq\.35[[:space:]]*$/a \\n## EUR\nN_no_favorable_lifestyle.category = sum(dat_all$admixture == "EUR", na.rm = TRUE)\naf_by_no_favorable_lifestyle.category.EUR = (N_all.EUR - N_no_favorable_lifestyle.category) / N_all.EUR\naf_by_no_favorable_lifestyle.category.EUR <- round(af_by_no_favorable_lifestyle.category.EUR,3)\naf_by_no_favorable_lifestyle.category.EUR\n\n## AFR\nN_no_favorable_lifestyle.category = sum(dat_all$admixture == "AFR", na.rm = TRUE)\naf_by_no_favorable_lifestyle.category.AFR = (N_all.AFR - N_no_favorable_lifestyle.category) / N_all.AFR\naf_by_no_favorable_lifestyle.category.AFR <- round(af_by_no_favorable_lifestyle.category.AFR,3)\naf_by_no_favorable_lifestyle.category.AFR' {} +


## add EUR and AFR to final result line
find . -type f -name '*model_fit*' -exec sed -i 's/age\.gteq = c(af_by_rt\.gteq\.35, af_by_tx\.gteq\.35, af_by_tx\.rt\.gteq\.35, af_by_prs\.gteq\.35, af_by_no_favorable_lifestyle\.category\.gteq\.35, af_by_combined\.gteq\.35)/age.gteq = c(af_by_rt.gteq.35, af_by_tx.gteq.35, af_by_tx.rt.gteq.35, af_by_prs.gteq.35, af_by_no_favorable_lifestyle.category.gteq.35, af_by_combined.gteq.35),\n  EUR = c(af_by_rt.EUR, af_by_tx.EUR, af_by_tx.rt.EUR, af_by_prs.EUR, af_by_no_favorable_lifestyle.category.EUR, af_by_combined.EUR),\n  AFR = c(af_by_rt.AFR, af_by_tx.AFR, af_by_tx.rt.AFR, af_by_prs.AFR, af_by_no_favorable_lifestyle.category.AFR, af_by_combined.AFR)/' {} +


# ## Replace admixture files in CCSS
# ## Admixture classification
# admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
# admixture$INDIVIDUAL <- sapply(strsplit(admixture$INDIVIDUAL,"_"), `[`, 1)



## Replace V18b with V19b
# find . -type f -name '*model_fit*' -exec sed -i 's/V18b/V19b/g' {} +
find . -type f -exec sed -i 's/V18b/V19b/g' {} +


## replace $sjlid with $ccssid
/home/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/V19b/19b/CCSS
find . -type f -name '*model_fit*' -exec sed -i 's/\$sjlid/\$ccssid/g' {} +


## uncomment 
find . -type f -name '*model_fit*' -exec sed -i 's/# EAS + AFR +/EAS + AFR +/' {} +


find . -type f -name '*model_fit*' -exec sed -i 's|admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_4481_Admixture_PCA_ethnicity.csv", sep = "\\t", header = T)|admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)|' {} +




## Replace in 20b
find . -type f -exec sed -i 's/V19b/V20b/g' {} +



find . -type f -name '*model_fit*' -exec sed -i 's/tertile\.Rdata/lt60.Rdata/g' {} +
find . -type f -name '*model_fit*' -exec sed -i 's/HEI2015_TOTAL_SCORE\.tertile\.category/HEI2015_TOTAL_SCORE.lt60.category/g' {} +
find . -type f -name '*model_fit*' -exec sed -i 's/category = "3rd"/category = "No"/g' {} +


## Replace Current smoker with Smoker ever variable
cd /home/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/V21b/21b
find . -type f  -exec sed -i 's/Current_smoker_yn/Smoker_ever_yn/g' {} +
find . -type f  -exec sed -i 's/smoker_former_or_never_yn_agesurvey/Smoker_ever_yn_agesurvey/g' {} +
find . -type f -exec sed -i 's/V20b/V21b/g' {} +
