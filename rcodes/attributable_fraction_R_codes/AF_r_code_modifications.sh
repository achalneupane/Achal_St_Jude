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