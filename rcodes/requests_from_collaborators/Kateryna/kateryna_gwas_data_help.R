library(haven)
library(dplyr)
# 1. in the ctcaegrades dataset, get the records with condition == Dyslipidemia - Hypercholerolemia or Dyslipidemia - Hypertriglyceridemia, MRN in the GWAS data.
# 2. For each MRN, only keep the records with gradedt <= labdt + 180. gradedt is in the ctcaegrades data and labdt is the lipid panel date in the gwas data.

# 3. For each mrn, if max(grade)>=2, then Dyslipidemia = 1


# CTCAE <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
gwas_data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/Dyslipidemia/GWAS/gwas_data.csv", header = T, sep = ",")

sum(grepl("Dyslipidemia - Hypercholerolemia|Dyslipidemia - Hypertriglyceridemia", CTCAE$condition))
# 11203


## ectract Dyslipidemia - Hypercholerolemia|Dyslipidemia - Hypertriglyceridemia from CTCAE
CTCAE.dyslipedemia <- CTCAE[grepl("Dyslipidemia - Hypercholesterolemia|Dyslipidemia - Hypertriglyceridemia", CTCAE$condition),]

# get the labdt to CTCAE extracted
CTCAE.dyslipedemia$labdt <- gwas_data$labdt[match(CTCAE.dyslipedemia$MRN, gwas_data$mrn)]


CTCAE.dyslipedemia$gradedt <- as.Date(CTCAE.dyslipedemia$gradedt, format = "%Y-%m-%d")
CTCAE.dyslipedemia$labdt <- as.Date(CTCAE.dyslipedemia$labdt, format = "%Y-%m-%d")


# For each MRN, only keep the records with gradedt <= labdt + 180
CTCAE.dyslipedemia_filtered <- CTCAE.dyslipedemia[CTCAE.dyslipedemia$gradedt <= CTCAE.dyslipedemia$labdt + 180, ]

# grade -9 = NA
# CTCAE.dyslipedemia_filtered$grade[CTCAE.dyslipedemia_filtered$grade == "-9"] <- NA

# samples with no MRN removed. not part of gwas_data
CTCAE.dyslipedemia_filtered <- CTCAE.dyslipedemia_filtered[!is.na(CTCAE.dyslipedemia_filtered$MRN),]



# 3. For each mrn, if max(grade)>=2, then Dyslipidemia = 1
CTCAE.dyslipedemia_filtered <- CTCAE.dyslipedemia_filtered %>%
  dplyr::group_by(MRN) %>%
  dplyr::mutate(Dyslipidemia = ifelse(all(is.na(grade)), NA, ifelse(max(as.numeric(grade), na.rm = TRUE) >= 2, 1, 0))) %>%
  dplyr::ungroup()

gwas_data$labdt <- as.Date(gwas_data$labdt)
# i=1

gwas_data$max_grade <- NA
gwas_data$max_grade_condition <- NA
gwas_data$sjlid <- NA
for (i in 1:nrow(gwas_data)) {
  cc <- subset(CTCAE.dyslipedemia_filtered, gwas_data$mrn[i] == CTCAE.dyslipedemia_filtered$MRN)
  cc.saved <- cc
  cc <- cc[which((abs(gwas_data$labdt[i] - cc$gradedt) <= 180/365.25)),]
  
max_grade_index <- which.max(cc$grade)[1]
if (!is.na(max_grade_index)) {
      gwas_data$max_grade[i] <- cc$grade[max_grade_index]
      gwas_data$max_grade_condition[i] <- cc$condition[max_grade_index]
      gwas_data$sjlid[i] <- cc$sjlid[1]
    } else {
      gwas_data$max_grade[i] <- "-9"
      gwas_data$max_grade_condition[i] <- NA
      gwas_data$sjlid[i] <- cc.saved$sjlid[1]
    }
}

gwas_data$status <- gwas_data$max_grade
gwas_data$status[gwas_data$max_grade>=2] <- 1
gwas_data$status[ gwas_data$max_grade>=0 & gwas_data$max_grade <2] <- 0

write.table(gwas_data, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/Dyslipidemia/GWAS/gwas_data_AN.txt", sep = "\t", col.names = T, row.names = F, quote = F)
