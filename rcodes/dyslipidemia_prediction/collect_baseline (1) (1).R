rm(list = ls())
library(haven)
library(dplyr)
library(janitor)
library(lubridate)
library(stringi)
sjlife_path <- "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/"

## earliest body measurement date is 2007 -- not available!!

function_combo_basic <- read_sas(paste0(sjlife_path, "Clinical data/function_combo_basic.sas7bdat")) %>%
  rename_with(tolower) %>%
  select(mrn, assmntdate, bmi, bmiadj) %>%  ## what is bmiadj?
  mutate(bmi = as.numeric(bmi), bmiadj = as.numeric(bmiadj)) %>% 
  filter(!is.na(bmi))

chemosum_dose <- read_sas(paste0(sjlife_path, "Clinical data/chemosum_dose.sas7bdat")) %>%
  rename_with(tolower) %>%
  select(mrn, anthracyclines_dose_5, carboplatin_dose_5, 
         cisplatin_dose_5) 

chemosum_yn <- read_sas(paste0(sjlife_path, "Clinical data/chemosum_yn.sas7bdat")) %>%
  rename_with(tolower) %>%
  select(mrn, corticosteroids_5) 

radiationsum_yn <- read_sas(paste0(sjlife_path, "Clinical data/radiationsum_yn.sas7bdat")) %>%
  rename_with(tolower) %>%
  select(mrn, cranial_5, anyrt_5) 

## question: variables to use about pelvic and abominal RT
radiation_dosimetry <- read_sas(paste0(sjlife_path, "Clinical data/radiation_dosimetry.sas7bdat")) %>%
  rename_with(tolower) %>%
  select(mrn, tbi, pelvis, abdomen, pancavgdose, maxabdrtdose, tbidose) %>%
  mutate(pancavgdose = ifelse(pancavgdose >= 10000, NA, pancavgdose),
         maxabdrtdose = ifelse(maxabdrtdose >= 10000, NA, maxabdrtdose),
         tbidose = ifelse(tbidose >= 10000, NA, tbidose)) %>%
  mutate(tbi = ifelse(tbi == "Yes", 1,
                      ifelse(tbi == "No", 0, NA)),
         pelvis = ifelse(pelvis == "Yes", 1,
                         ifelse(pelvis == "No", 0, NA)),
         abdomen = ifelse(abdomen == "Yes", 1,
                          ifelse(abdomen == "No", 0, NA)))


adult_healthhabits <- read_sas(paste0(sjlife_path, "Survey data/adult_healthhabits.sas7bdat")) %>%
  rename_with(tolower) %>%
  select(mrn, datecomp, evsm, cigmo, smnow, nopa,
         vpa10, vpadays, vpamin, mpa10, mpadays, mpamin,
         lpa10, lpadays, lpamin, ltpaw, wtlt, yoga, pa20) %>%
  filter(!is.na(datecomp))

adolescent_healthhabits <- read_sas(paste0(sjlife_path, "Survey data/adolescent_healthhabits.sas7bdat")) %>%
  rename_with(tolower) %>%
  select(mrn, datecomp, relation, vpa10, vpadays, vpamin, 
         mpa10, mpadays, mpamin,
         lpa10, lpadays, lpamin,
         activity, play60, pa20days, patone) %>%
  filter(!is.na(datecomp) & relation == 2) %>%
  arrange(mrn, datecomp, vpa10, vpadays, vpamin, mpa10, mpadays, mpamin, lpa10,
          lpadays, lpamin, activity, play60, pa20days, patone) %>%
  group_by(mrn, datecomp) %>%
  filter(row_number() == 1)

ctcaegrades <- read_sas(paste0(sjlife_path, "Event data/ctcaegrades.sas7bdat"))%>%
  rename_with(tolower) %>%
  select(mrn, condition, gradedt, grade, category) %>%
  filter(grade != -9, !is.na(gradedt)) ## mrn 2224 has one missing grade date for diabetes, removed

ctcaegrades[360047, "gradedt"] <- as.Date("2012-06-04") ## correct one error in the grade date

dysl_pop <- readRDS("dysl_pop_sa.rds")


## chemotherapy doses
dysl_chemo <- dysl_pop %>%
  select(mrn) %>%
  left_join(chemosum_dose, by = "mrn") %>%
  left_join(chemosum_yn, by = "mrn")



## prediction data = last contact date - 15 years
dysl_pop$predictdt <- dysl_pop$lstcondt %m-% years (15)
dysl_pop$agepredict <- with(dysl_pop,
                            time_length(predictdt - dob, "year"))

dysl_pop$adult <- as.numeric(dysl_pop$agepredict >= 18)

## at the age of prediction: 1065 adolescents, 1494 adults
## bmi_diff_dt: 1001 - 5479
dysl_bmi <- function_combo_basic %>% 
  right_join(select(dysl_pop, mrn, predictdt), by = "mrn") %>%
  mutate(diffdt = as.numeric(assmntdate - predictdt)) %>%
  group_by(mrn) %>%
  filter(diffdt == min(diffdt)) %>%
  summarise(mrn = first(mrn), 
            bmiadj = first(bmiadj))


## smoking status
## adolescents: do not consider
## adult: classified as current smoker, former smoker, never smoker or unknown
## issue: the earliest smoking status is 944 days after the prediction time
dysl_smoker <- dysl_pop %>% filter(adult == 1) %>% select(mrn) %>%
  left_join(adult_healthhabits %>% select(mrn,datecomp, evsm, cigmo, smnow), by = "mrn") %>%
  mutate(smoker = ifelse(is.na(evsm) & is.na(cigmo) & is.na(smnow),
                         NA, ifelse(cigmo == 2 & evsm == 1, "Former smoke",
                                    ifelse(cigmo == 1, "Now smoke",
                                           ifelse(cigmo == 2 & evsm == 2, "Never smoke",
                                                  NA))))) %>%
  left_join(dysl_pop %>% select(mrn, predictdt), by = "mrn") %>%
  mutate(diffdt = as.numeric(datecomp - predictdt)) %>%
  group_by(mrn) %>%
  arrange(diffdt) %>%
  summarise(mrn = first(mrn), smoker = first(smoker)) %>%
  right_join(dysl_pop %>% select(mrn), by = "mrn")


## compute MET for each adolescent with MRN in the study population
pa_adoles <- adolescent_healthhabits %>%
  select(mrn, datecomp, vpa10, vpadays, vpamin,
         mpa10, mpadays, mpamin, lpa10, lpadays, lpamin,
         activity, play60, pa20days, patone) %>%
  mutate(wvpa = NA, wmpa = NA, mvpawk = NA, wmspa = NA, mets = NA)


for (i in 1:nrow(pa_adoles)) {
  if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1 & is.na(pa_adoles$vpadays[i])) {
    pa_adoles$vpadays[i] <- 1
  }
  
  if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1 & is.na(pa_adoles$vpamin[i])) {
    pa_adoles$vpamin[i] <- 10
  }
  
  if (!is.na(pa_adoles$vpamin[i]) & pa_adoles$vpamin[i] > 360) {
    pa_adoles$vpamin[i] <- 360
  }
  
  if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1) {
    pa_adoles$wvpa[i] = pa_adoles$vpadays[i] * pa_adoles$vpamin[i]
  }
  
  if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2) {
    pa_adoles$wvpa[i] = 0
  }
  
  if (is.na(pa_adoles$wvpa[i])) {
    if ((is.na(pa_adoles$vpa10[i]) | (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
        (!is.na(pa_adoles$activity[i]) & pa_adoles$activity[i] == 2)) {
      pa_adoles$wvpa[i] <- 0
    } 
    
    if ((is.na(pa_adoles$vpa10[i])|(!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
        !(is.na(pa_adoles$play60[i])|(!is.na(pa_adoles$play60[i]) & pa_adoles$play60[i] == 1))) {
      pa_adoles$wvpa[i] <- 0
    }
    
    if ((is.na(pa_adoles$vpa10[i]) | (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
        !(is.na(pa_adoles$pa20days[i]) | (!is.na(pa_adoles$pa20days[i]) & pa_adoles$pa20days[i] == 1))) {
      pa_adoles$wvpa[i] <- (pa_adoles$pa20days[i] - 1) * 20
    }
  }
  
  if (is.na(pa_adoles$mpa10[i]) & 
      ((!is.na(pa_adoles$mpadays[i]))|(!is.na(pa_adoles$mpamin[i])))) {
    pa_adoles$mpa10[i] <- 1
  }
  
  if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1 & is.na(pa_adoles$mpadays[i])) {
    pa_adoles$mpadays[i] <- 1
  }
  
  if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1 & is.na(pa_adoles$mpamin[i])) {
    pa_adoles$mpamin[i] <- 10
  }
  
  if (!is.na(pa_adoles$mpamin[i]) & pa_adoles$mpamin[i] > 360) pa_adoles$mpamin[i] <- 360
  if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1) {
    pa_adoles$wmpa[i] <- 
      pa_adoles$mpadays[i] * pa_adoles$mpamin[i]
  }
  if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 2) {
    pa_adoles$wmpa[i] <-0
  }
  
  if (is.na(pa_adoles$wmpa[i])) {
    if ((is.na(pa_adoles$mpa10[i]) | pa_adoles$mpa10[i] == 2) &
        (!is.na(pa_adoles$activity[i]) & pa_adoles$activity[i] == 2)) {
      pa_adoles$wmpa[i] <- 0
    } 
    
    if ((is.na(pa_adoles$mpa10[i]) | pa_adoles$mpa10[i] == 2) &
        !(is.na(pa_adoles$play60[i]) | pa_adoles$play60[i] == 1)) {
      pa_adoles$wmpa[i] <- (pa_adoles$play60[i] - 1) * 60
    } 
    
    if (!is.na(pa_adoles$wvpa[i])) pa_adoles$wmpa[i] <- 0
  }
  
  
  if (is.na(pa_adoles$wmpa[i]) & is.na(pa_adoles$wvpa[i])) {
    if (!is.na(pa_adoles$play60[i]) & pa_adoles$play60[i] == 1) {
      pa_adoles$wvpa[i] <- pa_adoles$wmpa[i] <- 0
    }
  }
  
  if (!is.na(pa_adoles$patone[i]) & pa_adoles$patone[i] > 1) {
    pa_adoles$wmspa[i] <- (pa_adoles$patone[i] - 1) * 20
  } else {
    pa_adoles$wmspa[i] <- 0
  }
  
  pa_adoles$mvpawk[i] <- pa_adoles$wmpa[i] + 2 * pa_adoles$wvpa[i] +
    2 * pa_adoles$wmspa[i]
  
  if (is.na(pa_adoles$wvpa[i])) pa_adoles$mvpawk[i] <- pa_adoles$wmpa[i]
  
  if (!is.na(pa_adoles$mvpawk[i]) & pa_adoles$mvpawk[i] > 2520) pa_adoles$mvpawk[i] <- 2520
  
  pa_adoles$mets[i] <- 3 * pa_adoles$mvpawk[i]
  
}

pa_adoles_final <- select(pa_adoles, mrn, datecomp, mets)

## compute adult mets
pa_adults <- adult_healthhabits %>%
  select(mrn, datecomp, nopa,
         vpa10, vpadays, vpamin, mpa10, mpadays, mpamin,
         lpa10, lpadays, lpamin, ltpaw, wtlt, yoga, pa20) %>%
  mutate(wvpa = NA, wmpa = NA, mvpawk = NA, mets = NA)


for (i in 1:nrow(pa_adults)) {
  if (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 1 & is.na(pa_adults$vpadays[i])) {
    pa_adults$vpadays[i] <- 1
  }
  
  if (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 1 & is.na(pa_adults$vpamin[i])) {
    pa_adults$vpamin[i] <- 10
  }
  
  if (!is.na(pa_adults$vpamin[i]) & pa_adults$vpamin[i] > 360) {
    pa_adults$vpamin[i] <- 360
  }
  
  if (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 1) {
    pa_adults$wvpa[i] <- pa_adults$vpadays[i] * pa_adults$vpamin[i]
  }
  
  if (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 2) {
    pa_adults$wvpa[i] <- 0
  }
  
  if (is.na(pa_adults$wvpa[i])) {
    if ((is.na(pa_adults$vpa10[i]) | (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 2)) &
        (!is.na(pa_adults$nopa[i]) & pa_adults$nopa[i] == 1) &
        !(is.na(pa_adults$pa20[i]) | (!is.na(pa_adults$pa20[i]) & pa_adults$pa20[i] == 0))) {
      pa_adults$wvpa[i] <- 0
    } 
    
    if ((is.na(pa_adults$vpa10[i])|(!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 2)) &
        !(is.na(pa_adults$pa20[i])|(!is.na(pa_adults$pa20[i]) & pa_adults$pa20[i] == 0))) {
      pa_adults$wvpa[i] <- pa_adults$pa20[i] * 20
    }
    
    if ((is.na(pa_adults$vpa10[i]) | (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 2)) &
        (!is.na(pa_adults$nopa[i]) & pa_adults$nopa[i] == 2) &
        (is.na(pa_adults$pa20[i]) | (!is.na(pa_adults$pa20[i]) & pa_adults$pa20[i] == 0))) {
      pa_adults$wvpa[i] <- 0
    }
  }
  
  if (is.na(pa_adults$mpa10[i]) & 
      ((!is.na(pa_adults$mpadays[i]))|(!is.na(pa_adults$mpamin[i])))) {
    pa_adults$mpa10[i] <- 1
  }
  
  if (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 1 & is.na(pa_adults$mpadays[i])) {
    pa_adults$mpadays[i] <- 1
  }
  
  if (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 1 & is.na(pa_adults$mpamin[i])) {
    pa_adults$mpamin[i] <- 10
  }
  
  if (!is.na(pa_adults$mpamin[i]) & pa_adults$mpamin[i] > 360) pa_adults$mpamin[i] <- 360
  if (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 1) {
    pa_adults$wmpa[i] <- 
      pa_adults$mpadays[i] * pa_adults$mpamin[i]
  }
  if (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 2) {
    pa_adults$wmpa[i] <-0
  }
  
  if (is.na(pa_adults$wmpa[i])) {
    if ((is.na(pa_adults$mpa10[i]) | 
         (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 2)) &
        (!is.na(pa_adults$nopa[i]) & pa_adults$nopa[i] == 1)) {
      pa_adults$wmpa[i] <- 0
    } 
    
    if (!is.na(pa_adults$wvpa[i])) {
      pa_adults$wmpa[i] <- 0
    } 
    
  }
  
  
  
  
  pa_adults$mvpawk[i] <- pa_adults$wmpa[i] + 2 * pa_adults$wvpa[i] 
  
  if (!is.na(pa_adults$mvpawk[i]) & pa_adults$mvpawk[i] > 2520) pa_adults$mvpawk[i] <- 2520
  
  pa_adults$mets[i] <- 3 * pa_adults$mvpawk[i]
  
}

pa_adults_final <- select(pa_adults, mrn, datecomp, mets)

## 944 to 5483
dysl_pa <- rbind(pa_adoles_final, pa_adults_final) %>%
  right_join(dysl_pop %>% select(mrn, predictdt), by = "mrn") %>%
  mutate(diffdt = as.numeric(datecomp - predictdt)) %>%
  arrange(mrn, diffdt) %>%
  group_by(mrn) %>%
  filter(row_number() == 1) %>%
  select(mrn, mets)




## radiation therapy info
dysl_rt <- dysl_pop %>%
  select(mrn) %>%
  left_join(radiationsum_yn, by = "mrn") %>%
  left_join(radiation_dosimetry, by = "mrn")


## comorbidity condition
## e.g. diabetes:
## diabetes = 1 and missing_diabetes = 0 if:
##    maximum of baseline grades >= 2;
## diabetes = 0.5 and missing_diabetes = 0 (unknown at baseline) if:
##    no baseline grades and earliest grade >= 2
## diabetes = 0 and missing_diabetes = 0 if:
##    either maximum of baseline grades == 0,1 or 
##    there are no baseline grades abd first diabetes grade == 0,1
## diabetes = 0 and missing_diabetes = 1 if:
##    no diabetes grades available for the subject



## following the cardiomyopathy paper, defined as CTCAE grade >= 2

## diabetes
diabetes_grades <- ctcaegrades %>% 
  filter(condition == "Abnormal glucose metabolism") %>%
  filter(mrn %in% dysl_pop$mrn) %>%
  left_join(dysl_pop[, c("mrn", "predictdt", "lstcondt")], by = "mrn")


base_diabetes <- diabetes_grades %>% 
  filter(gradedt <= predictdt) %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>% 
  summarise(baselinedt = last(gradedt), baselinegrade = max(grade))


base_diabetes_case <- base_diabetes %>% 
  filter(baselinegrade >= 2) %>% 
  select(mrn) %>%
  mutate(diabetes = 1, missing_diabetes = 0)

base_diabetes_ctrl <- base_diabetes %>%
  filter(baselinegrade < 2) %>%
  select(mrn) %>%
  mutate(diabetes = 0, missing_diabetes = 0)

earliest_diabetes <- diabetes_grades %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>%
  filter(gradedt == min(gradedt)) %>%
  rename(earliestdt = gradedt, earliestgrade = grade)

earliest_diabetes_ctrl <- earliest_diabetes %>%
  filter(earliestgrade <= 1 & !(mrn %in% base_diabetes$mrn)) %>%
  select(mrn) %>%
  mutate(diabetes = 0, missing_diabetes = 0)

diabetes_ukbase <- earliest_diabetes %>%
  filter(earliestgrade >= 2 & !(mrn %in% base_diabetes$mrn)) %>%
  select(mrn) %>%
  mutate(diabetes = 0.5, missing_diabetes = 0)


diabetes_rslt <- bind_rows(base_diabetes_case,
                           base_diabetes_ctrl,
                           earliest_diabetes_ctrl,
                           diabetes_ukbase) %>%
  right_join(dysl_pop %>% select(mrn), by = "mrn") %>%
  mutate(missing_diabetes = ifelse(is.na(missing_diabetes), 1, missing_diabetes),
         diabetes = ifelse(is.na(diabetes), 0, diabetes))

## growth hormone deficiency: grade >= 1
ghd_grades <- ctcaegrades %>% 
  filter(condition %in% c("Adult GHD", "Childhood GHD")) %>%
  filter(mrn %in% dysl_pop$mrn) %>%
  left_join(dysl_pop[, c("mrn", "predictdt", "lstcondt")], by = "mrn")


base_ghd <- ghd_grades %>% 
  filter(gradedt <= predictdt) %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>% 
  summarise(baselinedt = last(gradedt), baselinegrade = max(grade))


base_ghd_case <- base_ghd %>% 
  filter(baselinegrade >= 1) %>% 
  select(mrn) %>%
  mutate(ghd = 1, missing_ghd = 0)

base_ghd_ctrl <- base_ghd %>%
  filter(baselinegrade < 1) %>%
  select(mrn) %>%
  mutate(ghd = 0, missing_ghd = 0)

earliest_ghd <- ghd_grades %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>%
  filter(gradedt == min(gradedt)) %>%
  rename(earliestdt = gradedt, earliestgrade = grade)

earliest_ghd_ctrl <- earliest_ghd %>%
  filter(earliestgrade == 0 & !(mrn %in% base_ghd$mrn)) %>%
  select(mrn) %>%
  mutate(ghd = 0, missing_ghd = 0)

ghd_ukbase <- earliest_ghd %>%
  filter(earliestgrade >= 1 & !(mrn %in% base_ghd$mrn)) %>%
  select(mrn) %>%
  mutate(ghd = 0.5, missing_ghd = 0)


ghd_rslt <- bind_rows(base_ghd_case,
                      base_ghd_ctrl,
                      earliest_ghd_ctrl,
                      ghd_ukbase) %>%
  right_join(dysl_pop %>% select(mrn), by = "mrn") %>%
  mutate(missing_ghd = ifelse(is.na(missing_ghd), 1, missing_ghd),
         ghd = ifelse(is.na(ghd), 0, ghd))

## hypogonadism
hypogonadism_grades <- ctcaegrades %>% 
  filter(condition %in% c("Hypogonadism Central (LH/FSHD)")) %>%
  filter(mrn %in% dysl_pop$mrn) %>%
  left_join(dysl_pop[, c("mrn", "predictdt", "lstcondt")], by = "mrn")


base_hypogonadism <- hypogonadism_grades %>% 
  filter(gradedt <= predictdt) %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>% 
  summarise(baselinedt = last(gradedt), baselinegrade = max(grade))


base_hypogonadism_case <- base_hypogonadism %>% 
  filter(baselinegrade >= 2) %>% 
  select(mrn) %>%
  mutate(hypogonadism = 1, missing_hypogonadism = 0)

base_hypogonadism_ctrl <- base_hypogonadism %>%
  filter(baselinegrade < 2) %>%
  select(mrn) %>%
  mutate(hypogonadism = 0, missing_hypogonadism = 0)

earliest_hypogonadism <- hypogonadism_grades %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>%
  filter(gradedt == min(gradedt)) %>%
  rename(earliestdt = gradedt, earliestgrade = grade)

earliest_hypogonadism_ctrl <- earliest_hypogonadism %>%
  filter(earliestgrade < 2 & !(mrn %in% base_hypogonadism$mrn)) %>%
  select(mrn) %>%
  mutate(hypogonadism = 0, missing_hypogonadism = 0)

hypogonadism_ukbase <- earliest_hypogonadism %>%
  filter(earliestgrade >= 2 & !(mrn %in% base_hypogonadism$mrn)) %>%
  select(mrn) %>%
  mutate(hypogonadism = 0.5, missing_hypogonadism = 0)


hypogonadism_rslt <- bind_rows(base_hypogonadism_case,
                               base_hypogonadism_ctrl,
                               earliest_hypogonadism_ctrl,
                               hypogonadism_ukbase) %>%
  right_join(dysl_pop %>% select(mrn), by = "mrn") %>%
  mutate(missing_hypogonadism = ifelse(is.na(missing_hypogonadism), 1, missing_hypogonadism),
         hypogonadism = ifelse(is.na(hypogonadism), 0, hypogonadism))

## hypothyroidism
hyperthyroidism_grades <- ctcaegrades %>% 
  filter(condition %in% c("Hypothyroidism Central", "Hypothyroidism Primary")) %>%
  filter(mrn %in% dysl_pop$mrn) %>%
  left_join(dysl_pop[, c("mrn", "predictdt", "lstcondt")], by = "mrn")


base_hyperthyroidism <- hyperthyroidism_grades %>% 
  filter(gradedt <= predictdt) %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>% 
  summarise(baselinedt = last(gradedt), baselinegrade = max(grade))


base_hyperthyroidism_case <- base_hyperthyroidism %>% 
  filter(baselinegrade >= 2) %>% 
  select(mrn) %>%
  mutate(hyperthyroidism = 1, missing_hyperthyroidism = 0)

base_hyperthyroidism_ctrl <- base_hyperthyroidism %>%
  filter(baselinegrade < 2) %>%
  select(mrn) %>%
  mutate(hyperthyroidism = 0, missing_hyperthyroidism = 0)

earliest_hyperthyroidism <- hyperthyroidism_grades %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>%
  filter(gradedt == min(gradedt)) %>%
  summarise(mrn = first(mrn), earliestdt = first(gradedt), 
            earliestgrade = max(grade))


earliest_hyperthyroidism_ctrl <- earliest_hyperthyroidism %>%
  filter(earliestgrade < 2 & !(mrn %in% base_hyperthyroidism$mrn)) %>%
  select(mrn) %>%
  mutate(hyperthyroidism = 0, missing_hyperthyroidism = 0)

hyperthyroidism_ukbase <- earliest_hyperthyroidism %>%
  filter(earliestgrade >= 2 & !(mrn %in% base_hyperthyroidism$mrn)) %>%
  select(mrn) %>%
  mutate(hyperthyroidism = 0.5, missing_hyperthyroidism = 0)


hyperthyroidism_rslt <- bind_rows(base_hyperthyroidism_case,
                                  base_hyperthyroidism_ctrl,
                                  earliest_hyperthyroidism_ctrl,
                                  hyperthyroidism_ukbase) %>%
  right_join(dysl_pop %>% select(mrn), by = "mrn") %>%
  mutate(missing_hyperthyroidism = ifelse(is.na(missing_hyperthyroidism), 1, missing_hyperthyroidism),
         hyperthyroidism = ifelse(is.na(hyperthyroidism), 0, hyperthyroidism))

## Cardiomyopathy
cardiomyopathy_grades <- ctcaegrades %>% 
  filter(condition %in% c("Cardiomyopathy")) %>%
  filter(mrn %in% dysl_pop$mrn) %>%
  left_join(dysl_pop[, c("mrn", "predictdt", "lstcondt")], by = "mrn")


base_cardiomyopathy <- cardiomyopathy_grades %>% 
  filter(gradedt <= predictdt) %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>% 
  summarise(baselinedt = last(gradedt), baselinegrade = max(grade))


base_cardiomyopathy_case <- base_cardiomyopathy %>% 
  filter(baselinegrade >= 2) %>% 
  select(mrn) %>%
  mutate(cardiomyopathy = 1, missing_cardiomyopathy = 0)

base_cardiomyopathy_ctrl <- base_cardiomyopathy %>%
  filter(baselinegrade < 2) %>%
  select(mrn) %>%
  mutate(cardiomyopathy = 0, missing_cardiomyopathy = 0)

earliest_cardiomyopathy <- cardiomyopathy_grades %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>%
  filter(gradedt == min(gradedt)) %>%
  rename(earliestdt = gradedt, earliestgrade = grade)

earliest_cardiomyopathy_ctrl <- earliest_cardiomyopathy %>%
  filter(earliestgrade < 2 & !(mrn %in% base_cardiomyopathy$mrn)) %>%
  select(mrn) %>%
  mutate(cardiomyopathy = 0, missing_cardiomyopathy = 0)

cardiomyopathy_ukbase <- earliest_cardiomyopathy %>%
  filter(earliestgrade >= 2 & !(mrn %in% base_cardiomyopathy$mrn)) %>%
  select(mrn) %>%
  mutate(cardiomyopathy = 0.5, missing_cardiomyopathy = 0)


cardiomyopathy_rslt <- bind_rows(base_cardiomyopathy_case,
                                 base_cardiomyopathy_ctrl,
                                 earliest_cardiomyopathy_ctrl,
                                 cardiomyopathy_ukbase) %>%
  right_join(dysl_pop %>% select(mrn), by = "mrn") %>%
  mutate(missing_cardiomyopathy = ifelse(is.na(missing_cardiomyopathy), 1, missing_cardiomyopathy),
         cardiomyopathy = ifelse(is.na(cardiomyopathy), 0, cardiomyopathy))


## Cardiomyopathy
ckd_grades <- ctcaegrades %>% 
  filter(condition %in% c("Chronic kidney disease")) %>%
  filter(mrn %in% dysl_pop$mrn) %>%
  left_join(dysl_pop[, c("mrn", "predictdt", "lstcondt")], by = "mrn")


base_ckd <- ckd_grades %>% 
  filter(gradedt <= predictdt) %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>% 
  summarise(baselinedt = last(gradedt), baselinegrade = max(grade))


base_ckd_case <- base_ckd %>% 
  filter(baselinegrade >= 1) %>% 
  select(mrn) %>%
  mutate(ckd = 1, missing_ckd = 0)

base_ckd_ctrl <- base_ckd %>%
  filter(baselinegrade < 1) %>%
  select(mrn) %>%
  mutate(ckd = 0, missing_ckd = 0)

earliest_ckd <- ckd_grades %>%
  select(mrn, gradedt, grade) %>%
  group_by(mrn) %>%
  filter(gradedt == min(gradedt)) %>%
  rename(earliestdt = gradedt, earliestgrade = grade)

earliest_ckd_ctrl <- earliest_ckd %>%
  filter(earliestgrade < 1 & !(mrn %in% base_ckd$mrn)) %>%
  select(mrn) %>%
  mutate(ckd = 0, missing_ckd = 0)

ckd_ukbase <- earliest_ckd %>%
  filter(earliestgrade >= 1 & !(mrn %in% base_ckd$mrn)) %>%
  select(mrn) %>%
  mutate(ckd = 0.5, missing_ckd = 0)


ckd_rslt <- bind_rows(base_ckd_case,
                      base_ckd_ctrl,
                      earliest_ckd_ctrl,
                      ckd_ukbase) %>%
  right_join(dysl_pop %>% select(mrn), by = "mrn") %>%
  mutate(missing_ckd = ifelse(is.na(missing_ckd), 1, missing_ckd),
         ckd = ifelse(is.na(ckd), 0, ckd))


final_data <- dysl_pop %>% select(mrn, dysl, agedx, gender, agepredict, adult) %>%
  left_join(dysl_bmi, by = "mrn") %>%
  left_join(dysl_smoker, by = "mrn") %>%
  left_join(dysl_pa, by = "mrn") %>%
  left_join(diabetes_rslt %>% select(mrn, diabetes), by = "mrn") %>%
  left_join(ghd_rslt %>% select(mrn, ghd), by = "mrn") %>%
  left_join(hypogonadism_rslt %>% select(mrn, hypogonadism), by = "mrn") %>%
  left_join(hyperthyroidism_rslt %>% select(mrn, hyperthyroidism), by = "mrn") %>%
  left_join(cardiomyopathy_rslt %>% select(mrn, cardiomyopathy), by = "mrn") %>%
  left_join(ckd_rslt %>% select(mrn, ckd), by = "mrn") %>%
  left_join(dysl_chemo, by = "mrn") %>%
  left_join(dysl_rt, by = "mrn")

saveRDS(final_data, file = "results/dysl_analytic_data_sa.rds")
