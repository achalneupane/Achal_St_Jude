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