rm(list = ls())
library(haven)
library(dplyr)
library(janitor)
library(lubridate)
source("utilities.R")
sjlife_path <- "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/"

diagnosis <- read_sas(paste0(sjlife_path, "Clinical data/diagnosis.sas7bdat")) %>% 
  rename_with(tolower) %>%
  arrange(mrn)


demographics <- read_sas(paste0(sjlife_path, "Clinical data/demographics.sas7bdat")) %>%
  rename_with(tolower) %>%
  arrange(mrn)  ## 10103 records

tracking <- read_sas(paste0(sjlife_path, "Tracking data/tracking.sas7bdat"))%>%
  rename_with(tolower) %>%
  arrange(mrn)  ## 11003 records

ctcaegrades <- read_sas(paste0(sjlife_path, "Event data/ctcaegrades.sas7bdat"))%>%
  rename_with(tolower) %>%
  select(mrn, condition, gradedt, grade, category) %>%
  filter(grade != -9) ## nonmissing grades

dyslipidemia_grades <- ctcaegrades %>%
  filter(condition %in% c("Dyslipidemia - Hypercholesterolemia",
                          "Dyslipidemia - Hypertriglyceridemia"))

cardiomyopathy_grades <- ctcaegrades %>%
  filter(condition %in% "Cardiomyopathy")

mi_grades <- ctcaegrades %>%
  filter(condition %in% "Myocardial Infarction")

cva_grades <- ctcaegrades %>%
  filter(condition %in% c("Cerebrovascular accident"))

pvd_grades <- ctcaegrades %>%
  filter(condition %in% "Vascular Disease")
## 14978 records on 4838 survivors
freeze_dt <- as.Date("2020-04-30")

lstcontact <- read_sas(paste0(sjlife_path, "Tracking data/lstcondt.sas7bdat"))%>%
  rename_with(tolower) %>%
  arrange(mrn) %>%
  mutate(lstcondt = ifelse(is.na(lstcondt), freeze_dt, lstcondt)) %>% ## if missing, set to be data freeze date
  mutate(lstcondt = as.Date(lstcondt))

## get patients' primary cancer diagnosis's visit info
diagnosis_prim <- diagnosis %>% filter(primdx == 1) ## 9232 records

## obtain information on those who finish the baseline visit
pop_1 <- left_join(tracking, 
                   select(demographics, -sjlid, -studypop, -abstraction, 
                          -abstractiontype, -abstraction_notes, -newstatus,
                          -deathdt),
                   by = "mrn") %>%
  left_join(select(diagnosis_prim,  -sjlid, -studypop, -abstraction, 
                   -abstractiontype, -abstraction_notes, -newstatus), 
            by = "mrn") %>%
  mutate(baselinevisitfinish = ifelse(is.na(baselinevisitfinish),
                                      freeze_dt,
                                      baselinevisitfinish),
         baselinevisitfinish = as.Date(baselinevisitfinish)) %>%
  mutate(age_base = (baselinevisitfinish - dob) / 365.25) %>%
  mutate(status = get_status(newstatus)) %>%
  mutate(age_base_2g = ifelse(age_base >= 18, 1, 0)) %>%
  mutate(studygrp = ifelse(studypop == "Survivor", 1, 2)) ## 10020 survivors, 983 community controls


## start with 10020 survivors 
pop_step1 <- filter(pop_1, studygrp == 1) ## 10020 selected


## exclude those status in (2, 3, 4, 5, 6)
pop_step2 <- filter(pop_step1, !status %in% c(2, 3, 4, 5, 6))  ## 6197 eligible survivors

## exclude nonparticipants
pop_step3 <- filter(pop_step2, !status %in% c(7, 8))  ## 5229 patients

## mrn of SJLIFE participants with baseline (even mild) dyslipidemia 
baseline_dyslipidemia <- dyslipidemia_grades %>% 
  inner_join(select(pop_step3, mrn, baselinevisitfinish), by = "mrn") %>%
  filter(grade >= 1) %>%
  filter(gradedt - baselinevisitfinish < 90)
  ## 1864 with baseline dyslipidemia
baseline_cardiomyopathy <- cardiomyopathy_grades %>% 
  inner_join(select(pop_step3, mrn, baselinevisitfinish), by = "mrn") %>%
  filter(grade >= 2) %>%
  filter(gradedt - baselinevisitfinish < 90)

baseline_mi <- mi_grades %>%
  inner_join(select(pop_step3, mrn, baselinevisitfinish), by = "mrn") %>%
  filter(grade >= 3) %>%
  filter(gradedt - baselinevisitfinish < 90)

baseline_cva <- cva_grades %>%
  inner_join(select(pop_step3, mrn, baselinevisitfinish), by = "mrn") %>%
  filter(grade >= 2) %>%
  filter(gradedt - baselinevisitfinish < 90)

baseline_pvd <- pvd_grades %>%
  inner_join(select(pop_step3, mrn, baselinevisitfinish), by = "mrn") %>%
  filter(grade >= 2) %>%
  filter(gradedt - baselinevisitfinish < 90)


pop_step4 <- pop_step3 %>%
  filter(!mrn %in% baseline_dyslipidemia$mrn,
         !mrn %in% baseline_cardiomyopathy$mrn,
         !mrn %in% baseline_cva$mrn,
         !mrn %in% baseline_pvd$mrn,
         !mrn %in% baseline_mi$mrn) ## 3010 participants

saveRDS(pop_step4, "results/dysl_pop.rds")

