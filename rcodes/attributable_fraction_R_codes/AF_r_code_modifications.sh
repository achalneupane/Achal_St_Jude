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
