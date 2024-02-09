library(haven)
library(dplyr)
# 1. in the ctcaegrades dataset, get the records with condition == Dyslipidemia - Hypercholerolemia or Dyslipidemia - Hypertriglyceridemia, MRN in the GWAS data.
# 2. For each MRN, only keep the records with gradedt <= labdt + 180. gradedt is in the ctcaegrades data and labdt is the lipid panel date in the gwas data.
# 3. For each mrn, if max(grade)>=2, then Dyslipidemia = 1


# CTCAE <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
gwas_data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/Dyslipidemia/GWAS/gwas_data.csv", header = T, sep = ",")

sum(grepl("Dyslipidemia - Hypercholerolemia|Dyslipidemia - Hypertriglyceridemia", CTCAE$condition))
# 11203

CTCAE.dyslipedemia <- CTCAE[grepl("Dyslipidemia - Hypercholerolemia|Dyslipidemia - Hypertriglyceridemia", CTCAE$condition),]

CTCAE.dyslipedemia$labdt <- gwas_data$labdt[match(CTCAE.dyslipedemia$MRN, gwas_data$mrn)]

CTCAE.dyslipedemia$gradedt <- as.Date(CTCAE.dyslipedemia$gradedt, format = "%Y-%m-%d")
CTCAE.dyslipedemia$labdt <- as.Date(CTCAE.dyslipedemia$labdt, format = "%Y-%m-%d")

CTCAE.dyslipedemia_filtered <- CTCAE.dyslipedemia[CTCAE.dyslipedemia$gradedt <= CTCAE.dyslipedemia$labdt + 180, ]

CTCAE.dyslipedemia_filtered$grade[CTCAE.dyslipedemia_filtered$grade == "-9"] <- NA

CTCAE.dyslipedemia_filtered <- CTCAE.dyslipedemia_filtered[!is.na(CTCAE.dyslipedemia_filtered$MRN),]

CTCAE.dyslipedemia_filtered <- CTCAE.dyslipedemia_filtered %>%
  dplyr::group_by(MRN) %>%
  dplyr::mutate(Dyslipidemia = ifelse(all(is.na(grade)), NA, ifelse(max(as.numeric(grade), na.rm = TRUE) >= 2, 1, 0))) %>%
  dplyr::ungroup()


cc <- cbind.data.frame(CTCAE.dyslipedemia_filtered$MRN, CTCAE.dyslipedemia_filtered$grade,CTCAE.dyslipedemia_filtered$gradedt, CTCAE.dyslipedemia_filtered$labdt, CTCAE.dyslipedemia_filtered$Dyslipidemia )

EUR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/rare_variants/SJLIFE_EURsamples.txt", header = T, sep = "\t")
AFR <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/rare_variants/SJLIFE_AFRsamples.txt", header = T, sep = "\t")                  

all.samples <- rbind.data.frame(EUR, AFR)

length(all.samples$IID)


