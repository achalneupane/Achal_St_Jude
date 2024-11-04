library(dplyr)
library(tidyr)
library(haven)

## Phenotype for all samples sequenced
all.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/sample_mapping.txt", header = F)

# CTCAE.data <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
sum(is.na(CTCAE.data$grade))
# 8075
# remove missing and unknown grades
CTCAE.data <- CTCAE.data[!is.na(CTCAE.data$grade),]
CTCAE.data <- CTCAE.data[-which(CTCAE.data$grade==-9),]

cc <- CTCAE.data[1:300, 1:11]

# CTCAE.data.1 <- CTCAE.data %>%
#   group_by(sjlid, condition, grade) %>%
#   arrange(desc(gradedt)) %>%
#   slice_head(n = 1) %>%
#   ungroup()



CTCAE.data.1 <- CTCAE.data %>%
  group_by(sjlid, condition) %>%
  mutate(all_zero = all(grade == 0),
         max_grade = max(grade),
         max_gradedt = min(gradedt[grade == max_grade])) %>%
  filter((all_zero & gradedt == max(gradedt)) | (grade == max_grade & gradedt == max_gradedt)) %>%
  arrange(sjlid, condition, desc(gradedt)) %>%
  select(-all_zero, -max_grade, -max_gradedt)


# # if all grade is zero, keep the first grade. If there ARE grades more than 0, then remove the rows with grade zero 
# CTCAE.data.2 <- CTCAE.data.1 %>%
#   group_by(sjlid, condition) %>%
#   mutate(all_zero = all(grade == 0)) %>%
#   filter(ifelse(all_zero, row_number() == 1, grade > 0)) %>%
#   arrange(sjlid, condition, desc(gradedt)) %>%
#   select(-all_zero)

# Keep the rows with maximum grade for each person and for each condition
CTCAE.data.2 <- CTCAE.data.1 %>%
  group_by(sjlid, condition) %>%
  mutate(all_zero = all(grade == 0),
         max_grade = ifelse(all_zero, max(grade), max(grade[grade > 0]))) %>%
  filter(grade == max_grade) %>%
  arrange(sjlid, condition, desc(gradedt)) %>%
  select(-all_zero, -max_grade)


all.conditions <- unique(CTCAE.data.2$condition)
i=1
all.conditions[i]

CTCAE.data.3 <- CTCAE.data.2[grepl(all.conditions[i], CTCAE.data.2$condition),]

# ## Now concat condition with grade
# CTCAE.data.2$condition_grade <- paste(CTCAE.data.2$condition, CTCAE.data.2$grade, sep = "_")
# CTCAE.data.2$condition_gradedt <- paste(CTCAE.data.2$condition, CTCAE.data.2$gradedt, sep = "_")
# CTCAE.data.2$condition_ageevent <- paste(CTCAE.data.2$condition, CTCAE.data.2$ageevent, sep = "_")

# cc <- CTCAE.data.2[1:10, 1:13]
# cc$condition_grade <- paste(cc$condition, cc$grade, sep = "_")
# cc$condition_gradedt <- paste(cc$condition, cc$gradedt, sep = "_")
# cc$condition_ageevent <- paste(cc$condition, cc$ageevent, sep = "_")

## Convert to wide format
CTCAE.data.4 <- CTCAE.data.3 %>%
  pivot_wider(
    id_cols = c(sjlid,studypop, sjlife_cohort, gender),
    names_from = c(condition, condition),
    values_from = c( grade, ageevent)
  )


CTCAE.data.3 <- cc %>%
  pivot_wider(
    id_cols = c(sjlid, studypop, sjlife_cohort, gender, sncount),
    names_from = c(condition_grade, condition_gradedt,condition_ageevent),
    values_from = c( grade, gradedt, ageevent)
  )

cc.2 <- CTCAE.data.1[CTCAE.data.1$sjlid == "SJL1527107",]
