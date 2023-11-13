rm(list = ls())
library(haven)
library(dplyr)
library(janitor)
library(lubridate)
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
  filter(condition %in% c("Dyslipidemia - Hypercholesterolemia",
                          "Dyslipidemia - Hypertriglyceridemia"))  %>%
  filter(grade != -9) ## nonmissing grades



## 14978 records on 4838 survivors
freeze_dt <- as.Date("2020-04-30")

lstcontact <- read_sas(paste0(sjlife_path, "Tracking data/lstcondt.sas7bdat"))%>%
  rename_with(tolower) %>%
  arrange(mrn) %>%
  mutate(lstcondt = ifelse(is.na(lstcondt), freeze_dt, lstcondt)) %>% ## if missing, set to be data freeze date
  mutate(lstcondt = as.Date(lstcondt))

## get patients' primary cancer diagnosis's visit info
diagnosis_prim <- diagnosis %>% filter(primdx == 1) ## 9232 records


get_status <- function(x) {
  y <- rep(NA, length(x))
  
  for (i in seq_along(x)) {
    if (x[i] == 3) {
      y[i] <- 1
    } else if (x[i] %in% c(8, 22)) {
      y[i] <- 2
    } else if (x[i] %in% c(15, 18)) {
      y[i] <- 3
    } else if (x[i] %in% c(1, 2, 13, 21, 11)) {
      y[i] <- 4
    } else if (x[i] %in% c(17, 20)) {
      y[i] <- 5
    } else if (x[i] %in% c(4, 5, 24)) {
      y[i] <- 6
    } else if (x[i] %in% c(9, 12, 19, 23)) {
      y[i] <- 7
    } else if (x[i] %in% c(6, 7, 10, 14, 99)) {
      y[i] <- 8
    }
  }
  return(y)
}

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


## get the last contact date, keep those with follow-up time longer than 15 years
pop_step3_lstcont <- left_join(pop_step3, 
                               lstcontact[, c("mrn", "lstcondt", "agelstcontact", "lstconsrc")], 
                               by = "mrn")


pop_step4 <- filter(pop_step3_lstcont,
                    time_length(lstcondt - diagdt, "year") >= 15)  ## 3937 subjects


## baseline dyslipidemia status (15 years before last condt)
baseline_dyslipidemia <- ctcaegrades %>% 
  inner_join(pop_step4, by = "mrn") %>%
  filter(gradedt <= lstcondt) %>%
  filter(time_length(lstcondt - gradedt, "year") >= 15) %>% ## 558 records on 278 survivors
  group_by(mrn) %>%
  summarise(base_dysl = as.numeric(max(grade) >= 2))  ## 131 with baseline dyslipidemia, 37 without

## having dyslipidemia grades within 15 years before last contact
dysl_grades <- ctcaegrades %>% 
  inner_join(pop_step4[, c("mrn", "lstcondt")], by = "mrn") %>%  ## 3916 patients with nonmissing dyslipidemia grades
  filter(gradedt <= lstcondt) %>%
  filter(time_length(lstcondt - gradedt, "year") <= 15) ## 16109records records on 3914 patients

dysl_rslt <- dysl_grades %>%
  group_by(mrn) %>%
  summarise(dysl = as.numeric(max(grade) >= 2))  ## 760 with dyslipidemia, 3174 without

pop_step5 <- left_join(dysl_rslt, baseline_dyslipidemia, by = "mrn") %>%
  mutate(base_dysl = ifelse(is.na(base_dysl), 0, base_dysl)) %>%
  filter(base_dysl == 0) ## 3783 participants, 633 with, 3150 without

pop_step5_final <- left_join(pop_step5[, c("mrn", "dysl")],
                             pop_step4[, c("mrn", "deathdt", "dob", "gender", "race",
                                           "ethnic", "racegrp", "racegrp2", "hispanic",
                                           "primdx", "diagdt", "agedx", "lstcondt", 
                                           "agelstcontact")],
                             by = "mrn")
saveRDS(pop_step5_final, "dysl_pop.rds")

## sensitivity analysis -- exclude baseline dyslipidemia grade = 1 and control with dyslipidemia grade = 1 at any time
base_dysl_sa <- ctcaegrades %>% 
  inner_join(pop_step4, by = "mrn") %>%
  filter(gradedt <= lstcondt) %>%
  filter(time_length(lstcondt - gradedt, "year") >= 15) %>% ## 334 records on 168 survivors
  group_by(mrn) %>%
  summarise(base_dysl = as.numeric(max(grade) >= 1))  ## 162 with baseline grade >= 1, 6 without

dysl_rslt_case <- filter(dysl_rslt, dysl == 1) ## case: max grade >= 2; 760 cases
dysl_rslt_ctrl <- dysl_grades %>%  ## 1934 controls
  group_by(mrn) %>%
  summarise(max_grade = max(grade)) %>%
  filter(max_grade == 0) %>%
  transmute(mrn = mrn, dysl = 0)

pop_step5_sa <- rbind(dysl_rslt_case, dysl_rslt_ctrl) %>%
  left_join(base_dysl_sa, by = "mrn") %>%
  mutate(base_dysl = ifelse(is.na(base_dysl), 0, base_dysl)) %>%
  filter(base_dysl == 0)  ## 627 cases; 1932 controls

pop_step5_sa_final <- left_join(pop_step5_sa[, c("mrn", "dysl")],
                                pop_step4[, c("mrn", "deathdt", "dob", "gender", "race",
                                              "ethnic", "racegrp", "racegrp2", "hispanic",
                                              "primdx", "diagdt", "agedx", "lstcondt", 
                                              "agelstcontact")],
                                by = "mrn")

saveRDS(pop_step5_sa_final, file = "dysl_pop_sa.rds")
