setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/nihms/")

####################
## Read geno file ##
####################
library(data.table)
library(haven)
raw <- as.data.frame(fread("subbed_recodeA.raw", header = T) )
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
colnames(raw) <- gsub("\\.", ":", HEADER)
sum(df$ID %in% colnames(raw))
# 67172

library(haven)
library(data.table)
CTCAE <- read_sas("/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
sum(grepl("cardio", CTCAE$organsys, ignore.case = T))
CTCAE <- CTCAE[grepl("cardio", CTCAE$organsys, ignore.case = T),]

## AA samples in SJLIFE
admixture <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/admixture.rds")
AFR <- admixture[admixture$ancestry =="AFR",]
AFR <- AFR$INDIVIDUAL[which(AFR$cohort == "Survivor")]

SJLIFE <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr18.preQC_biallelic_renamed_ID_edited.vcf.gz.fam", header = F)
SJLIFE.AFR <- SJLIFE[SJLIFE$V2 %in% AFR,]

## KEEP AFR only in CTCAE
CTCAE <- CTCAE[CTCAE$sjlid %in% AFR,]
table(CTCAE$condition)


conditions <- unique(CTCAE$condition)



library(dplyr)
library(tidyr)
CTCAE_transformed <- CTCAE %>%
  select(sjlid, condition) %>%  # Keep relevant columns
  distinct() %>%                # Remove duplicate rows
  mutate(value = 1) %>%         # Assign 1 for presence
  pivot_wider(names_from = condition, values_from = value, values_fill = list(value = 0))


unique(CTCAE$condition)



table(CTCAE$condition)



cc <- CTCAE[grepl("cardiomy", CTCAE$condition, ignore.case = T),]
cc <- cc[!duplicated(cc$sjlid),]
dim(cc)
