rm(list = ls())
library(haven)
library(dplyr)
library(janitor)
library(lubridate)
library(stringi)
library(readxl)
source("../utilities.R")

list_a <- read.table("data/ListA_APOB_LDLR_PCSK9_comb_carriers.txt",
                     header = T) %>%
  rename(sjlid = IID, list_a_carrier = Carrier_status)

list_ab <- read.table("data/ListAB_merged_comb_carriers.txt",
                      header = T) %>%
  rename(sjlid = IID, list_ab_carrier = Carrier_status)

list_abc <- read.table("data/ListABC_merged_comb_carriers.txt",
                       header = T) %>%
  rename(sjlid = IID, list_abc_carrier = Carrier_status)

sjlife_path <- "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/"

## earliest body measurement date is 2007 -- not available!!
tracking <- read_sas(paste0(sjlife_path, "Tracking data/tracking.sas7bdat"))%>%
  rename_with(tolower) %>%
  dplyr::select(mrn, sjlid)

demographics <- read_sas(paste0(sjlife_path, "Clinical data/demographics.sas7bdat")) %>%
  rename_with(tolower) %>%
  dplyr::select(mrn, dob, racegrp)

function_combo_basic <- read_sas(paste0(sjlife_path, "Clinical data/function_combo_basic.sas7bdat")) %>%
  rename_with(tolower) %>%
  dplyr::select(mrn, assmntdate, bmi, bmiadj) %>%  ## what is bmiadj?
  mutate(bmi = as.numeric(bmi), bmiadj = as.numeric(bmiadj)) %>% 
  filter(!is.na(bmi))

chemosum_dose <- read_sas(paste0(sjlife_path, "Clinical data/chemosum_dose.sas7bdat")) %>%
  rename_with(tolower) %>%
  dplyr::select(mrn, anthracyclines_dose_5, carboplatin_dose_5, 
                cisplatin_dose_5) 

chemosum_yn <- read_sas(paste0(sjlife_path, "Clinical data/chemosum_yn.sas7bdat")) %>%
  rename_with(tolower) %>%
  dplyr::select(mrn, corticosteroids_5) 

milli_labs <- read_sas(paste0(sjlife_path, "Clinical data/milli_labs.sas7bdat")) %>%
  rename_with(tolower) 

radiationsum_yn <- read_sas(paste0(sjlife_path, "Clinical data/radiationsum_yn.sas7bdat")) %>%
  rename_with(tolower) %>%
  dplyr::select(mrn, cranial_5, anyrt_5) 


radiation_dosimetry <- read_sas(paste0(sjlife_path, "Clinical data/radiation_dosimetry.sas7bdat")) %>%
  rename_with(tolower) %>%
  dplyr::select(mrn, maxabdrtdose, maxsegrtdose, brainrt_yn) %>%
  mutate(maxabdrtdose = ifelse(maxabdrtdose >= 10000, NA, maxabdrtdose))

adolescent_healthhabits <- read_sas(paste0(sjlife_path, "Survey data/adolescent_healthhabits.sas7bdat")) %>%
  rename_with(tolower)

adult_healthhabits <- read_sas(paste0(sjlife_path, "Survey data/adult_healthhabits.sas7bdat")) %>%
  rename_with(tolower)

adult_pa <- get_adult_pa() %>% filter(relation == 1) %>% ## for adult, use self-report
  dplyr::select(mrn, datecomp, mets)

adolescent_pa_all <- get_adolescent_pa() %>% ## drop self report, use parent report
  filter(relation %in% c(2, 3)) 

duplicated_adoles_pa_mrn <- adolescent_pa_all %>%
  group_by(mrn) %>%
  filter(n() >= 2)

adolescent_pa <- adolescent_pa_all %>%
  filter((!mrn %in% duplicated_adoles_pa_mrn) |
           (mrn %in% duplicated_adoles_pa_mrn &
              survey %in% c("11-17 Years of Age Parent Report",
                            "5-10 Years of Age Parent Report"))) %>%
  dplyr::select(mrn, datecomp, mets)


all_pa <- rbind(adult_pa, adolescent_pa) %>%
  arrange(mrn, datecomp) 

mets_mean <- mean(all_pa$mets, na.rm = T)
mets_sd <- sd(all_pa$mets, na.rm = T)

all_pa_z <- all_pa %>%
  mutate(mets_z = (mets - mets_mean) / mets_sd)


ctcaegrades <- read_sas(paste0(sjlife_path, "Event data/ctcaegrades.sas7bdat"))%>%
  rename_with(tolower) %>%
  dplyr::select(mrn, condition, gradedt, grade, category) %>%
  filter(grade != -9, !is.na(gradedt)) ## mrn 2224 has one missing grade date for diabetes, removed

ctcaegrades[360047, "gradedt"] <- as.Date("2012-06-04") ## correct one error in the grade date



## 14978 records on 4838 survivors
freeze_dt <- as.Date("2020-04-30")

lstcontact <- read_sas(paste0(sjlife_path, "Tracking data/lstcondt.sas7bdat"))%>%
  rename_with(tolower) %>%
  mutate(lstcondt = ifelse(is.na(lstcondt), freeze_dt, lstcondt)) %>% ## if missing, set to be data freeze date
  mutate(lstcondt = as.Date(lstcondt)) %>%
  transmute(mrn, lstcondt, death = as.numeric(lstconsrc %in% "Death Date"))


dysl_pop <- readRDS("results/dysl_pop.rds")


## chemotherapy doses
dysl_chemo <- dysl_pop %>%
  dplyr::select(mrn) %>%
  left_join(chemosum_dose, by = "mrn") %>%
  left_join(chemosum_yn, by = "mrn")



## bmi at baseline visit
dysl_bmi <- function_combo_basic %>% 
  right_join(dplyr::select(dysl_pop, mrn, baselinevisitfinish), by = "mrn") %>%
  filter(assmntdate < baselinevisitfinish + 180) %>%
  arrange(mrn, assmntdate) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, bmiadj)



## smoking status at baseline visit
## adolescents: counted as nonsmoker
## adult: classified as current smoker, former smoker, never smoker or unknown
## issue: the earliest smoking status is 944 days after the prediction time
smoker_dat <- adult_healthhabits %>%
  filter(relation == 1) %>%
  dplyr::select(mrn, datecomp, smoker) 

dysl_smoker <- match_visit(id = "mrn", 
                           data_a = dysl_pop %>% 
                             dplyr::select(mrn, baselinevisitfinish),
                           date_a = "baselinevisitfinish",
                           data_b = smoker_dat,
                           date_b = "datecomp") %>%
  left_join(demographics %>% dplyr::select(mrn, dob), by = "mrn") %>%
  mutate(baseage = yrdif(dob, baselinevisitfinish)) %>%
  mutate(smoker = ifelse(is.na(smoker), 4, smoker)) %>%
  dplyr::select(mrn, smoker)


## baseline PA
dysl_pa <- dysl_pop %>% dplyr::select(mrn, baselinevisitfinish) %>%
  left_join(all_pa_z, by = "mrn") %>%
  filter(datecomp < baselinevisitfinish + 180) %>%
  arrange(mrn, datecomp) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, mets, mets_z)




## radiation therapy info
dysl_rt <- dysl_pop %>%
  dplyr::select(mrn) %>%
  left_join(radiationsum_yn, by = "mrn") %>%
  left_join(radiation_dosimetry, by = "mrn") %>%
  mutate(brainrt_yn = ifelse(brainrt_yn == 1, 1, 0))

## comorbidity condition
## e.g. diabetes:
## diabetes = 1 and missing_diabetes = 0 if:
##    maximum of baseline grades >= 2;
## diabetes = 0.5 and missing_diabetes = 0 (unknown at baseline) if:
##    no baseline grades and earliest grade >= 2
## diabetes = 0 and missing_diabetes = 0 if:
##    either maximum of baseline grades == 0,1 or 
##    there are no baseline grades and first diabetes grade == 0,1
## diabetes = 0 and missing_diabetes = 1 if:
##    no diabetes grades available for the subject



## following the cardiomyopathy paper, defined as CTCAE grade >= 2

## diabetes
diabetes_grades <- ctcaegrades %>% 
  filter(condition == "Abnormal glucose metabolism") %>%
  filter(mrn %in% dysl_pop$mrn) 

base_diabetes <- dysl_pop %>%
  dplyr::select(mrn, baselinevisitfinish) %>%
  left_join(diabetes_grades, by = "mrn") %>%
  filter(gradedt < baselinevisitfinish + 180) %>%
  arrange(mrn, gradedt) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, grade) %>%
  right_join(dysl_pop %>% dplyr::select(mrn), by = "mrn") %>% 
  transmute(mrn = mrn, diabetes = as.numeric(grade %in% 2:5))



## Hypertension
htn_grades <- ctcaegrades %>% 
  filter(condition == "Hypertension") %>%
  filter(mrn %in% dysl_pop$mrn) 

base_htn <- dysl_pop %>%
  dplyr::select(mrn, baselinevisitfinish) %>%
  left_join(diabetes_grades, by = "mrn") %>%
  filter(gradedt < baselinevisitfinish + 180) %>%
  arrange(mrn, gradedt) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, grade) %>%
  right_join(dysl_pop %>% dplyr::select(mrn), by = "mrn") %>% 
  transmute(mrn = mrn, htn = as.numeric(grade %in% 2:5))

## growth hormone deficiency: grade >= 1
ghd_grades <- ctcaegrades %>% 
  filter(condition == "Hypertension") %>%
  filter(mrn %in% dysl_pop$mrn) 

base_ghd <- dysl_pop %>%
  dplyr::select(mrn, baselinevisitfinish) %>%
  left_join(diabetes_grades, by = "mrn") %>%
  filter(gradedt < baselinevisitfinish + 180) %>%
  arrange(mrn, gradedt) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, grade) %>%
  right_join(dysl_pop %>% dplyr::select(mrn), by = "mrn") %>% 
  transmute(mrn = mrn, ghd = as.numeric(grade %in% 1:5))

## hypogonadism
hypogonadism_grades <- ctcaegrades %>% 
  filter(condition == "Hypogonadism Central (LH/FSHD)") %>%
  filter(mrn %in% dysl_pop$mrn) 

base_hypogonadism <- dysl_pop %>%
  dplyr::select(mrn, baselinevisitfinish) %>%
  left_join(diabetes_grades, by = "mrn") %>%
  filter(gradedt < baselinevisitfinish + 180) %>%
  arrange(mrn, gradedt) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, grade) %>%
  right_join(dysl_pop %>% dplyr::select(mrn), by = "mrn") %>% 
  transmute(mrn = mrn, hypogonadism = as.numeric(grade %in% 2:5))

## hypothyroidism
hyperthyroidism_grades <- ctcaegrades %>% 
  filter(condition %in% c("Hypothyroidism Central", "Hypothyroidism Primary")) %>%
  filter(mrn %in% dysl_pop$mrn) 

base_hyperthyroidism <- dysl_pop %>%
  dplyr::select(mrn, baselinevisitfinish) %>%
  left_join(diabetes_grades, by = "mrn") %>%
  filter(gradedt < baselinevisitfinish + 180) %>%
  arrange(mrn, gradedt) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, grade) %>%
  right_join(dysl_pop %>% dplyr::select(mrn), by = "mrn") %>% 
  transmute(mrn = mrn, hyperthyroidism = as.numeric(grade %in% 2:5))

## Cardiomyopathy
cardiomyopathy_grades <- ctcaegrades %>% 
  filter(condition %in% c("Cardiomyopathy")) %>%
  filter(mrn %in% dysl_pop$mrn) 

base_cardiomyopathy <- dysl_pop %>%
  dplyr::select(mrn, baselinevisitfinish) %>%
  left_join(diabetes_grades, by = "mrn") %>%
  filter(gradedt < baselinevisitfinish + 180) %>%
  arrange(mrn, gradedt) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, grade) %>%
  right_join(dysl_pop %>% dplyr::select(mrn), by = "mrn") %>% 
  transmute(mrn = mrn, cardiomyopathy = as.numeric(grade %in% 2:5))

## Chronic kidney diseases
ckd_grades <- ctcaegrades %>% 
  filter(condition %in% c("Chronic kidney disease")) %>%
  filter(mrn %in% dysl_pop$mrn) 

base_ckd <- dysl_pop %>%
  dplyr::select(mrn, baselinevisitfinish) %>%
  left_join(diabetes_grades, by = "mrn") %>%
  filter(gradedt < baselinevisitfinish + 180) %>%
  arrange(mrn, gradedt) %>%
  group_by(mrn) %>%
  filter(row_number() == n()) %>%
  dplyr::select(mrn, grade) %>%
  right_join(dysl_pop %>% dplyr::select(mrn), by = "mrn") %>% 
  transmute(mrn = mrn, ckd = as.numeric(grade %in% 1:5))

## Outcome: dyslipidemia
dyslipidemia_grades <- ctcaegrades %>%
  filter(condition %in% c("Dyslipidemia - Hypercholesterolemia",
                          "Dyslipidemia - Hypertriglyceridemia")) %>%
  filter(grade >= 2) %>%
  mutate(dyslipidemia = 1) %>%
  dplyr::select(mrn, gradedt, dyslipidemia) %>%
  arrange(mrn, gradedt) %>%
  rename(dyslipidemia_dt = gradedt) %>%
  unique() %>%
  inner_join(dysl_pop %>% dplyr::select(mrn, baselinevisitfinish), 
             by = "mrn") %>%
  filter(dyslipidemia_dt > baselinevisitfinish) %>%
  dplyr::select(mrn, dyslipidemia, dyslipidemia_dt) %>%
  arrange(mrn, dyslipidemia_dt) %>%
  group_by(mrn) %>%
  filter(row_number() == 1)  

dysl_outcome <- dysl_pop %>% dplyr::select(mrn) %>%
  left_join(dyslipidemia_grades, by = "mrn") %>%
  mutate(dyslipidemia = ifelse(is.na(dyslipidemia), 0, dyslipidemia))

## Outcome: baseline lipids value
lipid_data <- milli_labs %>%
  dplyr::select(mrn, labdt, cholesterol_totl,
                hdl_cholesterol, ldl_cholesterol, triglyceride) %>%
  mutate(cholesterol_totl = as.numeric(cholesterol_totl),
         hdl_cholesterol = as.numeric(hdl_cholesterol),
         ldl_cholesterol = as.numeric(ldl_cholesterol),
         triglyceride = as.numeric(triglyceride)) %>%
  mutate(nonhdl_cholesterol = cholesterol_totl - hdl_cholesterol) %>%
  arrange(mrn, labdt) %>%
  group_by(mrn) %>%
  filter(row_number() == 1) 



ancestry <- read_xlsx("data/SJLIFE_genetic ancestry.xlsx") %>%
  rename_with(tolower) %>%
  dplyr::select(id, ancestry) %>%
  transmute(sjlid = id, white = as.numeric(ancestry == "EUR"))

prof677 <- read.table("data/PGS000677_GRCh38_prs38.profile", header = T) %>%
  rename_with(tolower) %>%
  dplyr::select(iid, score) %>%
  rename(sjlid = iid, score677 = score)

score677_mean <- mean(prof677$score677, na.rm = T)
score677_sd <- sd(prof677$score677, na.rm = T)

prof677_z <- mutate(prof677, score677_z = (score677 - score677_mean) / score677_sd)

prof686 <- read.table("data/PGS000686_GRCh38_prs38.profile", header = T) %>%
  rename_with(tolower) %>%
  dplyr::select(iid, score) %>%
  rename(sjlid = iid, score686 = score)


score686_mean <- mean(prof686$score686, na.rm = T)
score686_sd <- sd(prof686$score686, na.rm = T)

prof686_z <- mutate(prof686, score686_z = (score686 - score686_mean) / score686_sd)


prof688 <- read.table("data/PGS000688_GRCh38_prs38.profile", header = T) %>%
  rename_with(tolower) %>%
  dplyr::select(iid, score) %>%
  rename(sjlid = iid, score688 = score)

score688_mean <- mean(prof688$score688, na.rm = T)
score688_sd <- sd(prof688$score688, na.rm = T)

prof688_z <- mutate(prof688, score688_z = (score688 - score688_mean) / score688_sd)

prof699 <- read.table("data/PGS000699_GRCh38_prs38.profile", header = T) %>%
  rename_with(tolower) %>%
  dplyr::select(iid, score) %>%
  rename(sjlid = iid, score699 = score)

score699_mean <- mean(prof699$score699, na.rm = T)
score699_sd <- sd(prof699$score699, na.rm = T)

prof699_z <- mutate(prof699, score699_z = (score699 - score699_mean) / score699_sd)



prof888 <- read.table("data/PGS000888_GRCh38_prs38.profile", header = T) %>%
  rename_with(tolower) %>%
  dplyr::select(iid, score) %>%
  rename(sjlid = iid, score888 = score)

score888_mean <- mean(prof888$score888, na.rm = T)
score888_sd <- sd(prof888$score888, na.rm = T)

prof888_z <- mutate(prof888, score888_z = (score888 - score888_mean) / score888_sd)


final_data <- dysl_pop %>% 
  dplyr::select(mrn, baselinevisitfinish,dob, gender, racegrp, 
         agedx, diaggrp, deathdt, age_base) %>%
  left_join(lstcontact, by = "mrn") %>%
  left_join(dysl_bmi, by = "mrn") %>%
  left_join(dysl_smoker, by = "mrn") %>%
  left_join(dysl_pa, by = "mrn") %>%
  left_join(base_diabetes, by = "mrn") %>%
  left_join(base_ghd, by = "mrn") %>%
  left_join(base_hypogonadism, by = "mrn") %>%
  left_join(base_hyperthyroidism, by = "mrn") %>%
  left_join(base_cardiomyopathy, by = "mrn") %>%
  left_join(base_ckd, by = "mrn") %>%
  left_join(base_htn, by = "mrn") %>%
  left_join(dysl_chemo, by = "mrn") %>%
  left_join(dysl_rt, by = "mrn")  %>%
  left_join(dysl_outcome, by = "mrn") %>%
  left_join(lipid_data, by = "mrn") %>%
  left_join(tracking, by = "mrn") %>%
  left_join(ancestry, by = "sjlid") %>%
  left_join(prof677_z, by = "sjlid") %>%
  left_join(prof686_z, by = "sjlid") %>%
  left_join(prof688_z, by = "sjlid") %>%
  left_join(prof699_z, by = "sjlid") %>%
  left_join(prof888_z, by = "sjlid") %>%
  left_join(list_a, by = "sjlid") %>%
  left_join(list_ab, by = "sjlid") %>%
  left_join(list_abc, by = "sjlid")

saveRDS(final_data, file = "results/dysl_analytic_data.rds")

write.csv(final_data, file = "results/dyslipidemia_data.csv",
          row.names = F)
