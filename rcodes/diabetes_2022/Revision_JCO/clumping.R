library(dplyr)
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/Results_final_to_Cindy/Revision_task_JCO")

#########
## AFR ##
#########
afr.df <- read.table("TOP.AFR.only.with.P.5e-06.and.results.txt", header = T)
# colnames(afr.df)[1] <- "SNP"
lowest_p_variants <- afr.df %>%
  group_by(CHR) %>%
  slice(which.min(P))


# # index variants:
# # Index variants (EUR)
# chr19:22614854:A:G
# chr2:2434309:AAAATAAAATAAAATAAAATAAAAAC:A
# chr8:50638825:A:G
# # Index variants (AFR)
# chr5:9843018:T:C

index <- c("chr5:9843018:T:C")

df <- as.data.frame(afr.df$MarkerName)
colnames(df) <- "SNP"
# write.table(df, "AFR_5e-06_variants_clumping.txt", col.names = T, quote = F, row.names = F)

## Add AFR frequency
AFR.freq <- read.table("AFR.Alele1_freq_output_All.chr", header = T)
table(afr.df$MarkerName %in% AFR.freq$SNP )
afr.df$AFR.freq <- AFR.freq$MAF[match(afr.df$MarkerName, AFR.freq$SNP)]

## Add EUR frequency
EUR.freq <- read.table("AFR.Alele1_freq_output_All_in_replication_EUR.chr", header = T)
table(afr.df$MarkerName %in% EUR.freq$SNP )
afr.df$EUR.freq <- EUR.freq$MAF[match(afr.df$MarkerName, EUR.freq$SNP)]

## Add LD from Clumping
AFR.ld.clumping <- read.table("AFR_LD_clumped_variants_chr_ALL.ld", header = T)

afr.df$ld_with_index_var <- NULL
afr.df$index_var <- NULL
afr.df$index.var.pos <- NULL
for (chr in 1:22){
  print(chr)
  CHR= paste0("chr", chr, ":")
  index.var <- lowest_p_variants$MarkerName[grepl(CHR, lowest_p_variants$MarkerName)]
  index.var.pos <- lowest_p_variants$BP[grepl(CHR, lowest_p_variants$MarkerName)]
  afr.df$index_var[afr.df$CHR == chr] <- index.var
  afr.df$index.var.pos[afr.df$CHR == chr] <- index.var.pos
  AFR.ld.clumping.tmp <- AFR.ld.clumping[AFR.ld.clumping$CHR_A %in% chr,]
  if (nrow(AFR.ld.clumping.tmp) == 0) {
  next
  }
  
  # Swap SNP_A and SNP_B based on index.var to keep index always on the right column
  AFR.ld.clumping.tmp <- AFR.ld.clumping.tmp %>%
    mutate(
      temp_SNP = if_else(SNP_A == index.var, SNP_B, SNP_A),
      SNP_B = if_else(SNP_A == index.var, index.var, SNP_B),
      SNP_A = temp_SNP
    ) %>%
    select(-temp_SNP)
  AFR.ld.clumping.tmp <- AFR.ld.clumping.tmp[AFR.ld.clumping.tmp$SNP_B %in% index.var,]
  if (nrow(AFR.ld.clumping.tmp) == 0) {
    next
  }
  for (j in 1:nrow(AFR.ld.clumping.tmp)){
    other.var <- AFR.ld.clumping.tmp$SNP_A[j]
    afr.df$ld_with_index_var[afr.df$MarkerName == other.var] <- AFR.ld.clumping.tmp$R2[AFR.ld.clumping.tmp$SNP_A == other.var & AFR.ld.clumping.tmp$SNP_B == index.var]
}
}

afr.df$distance_from_index_BP <- afr.df$BP - afr.df$index.var.pos

write.table(afr.df, "AFR_clumping_results.txt", col.names = T, quote = F, row.names = F)

#########
## EUR ##
#########
eur.df <- read.table("TOP.EUR.only.with.P.5e-06.and.results.from.meta.variants.txt", header = T)

index <- c("chr2:2434309:AAAATAAAATAAAATAAAATAAAAAC:A", "chr8:50638825:A:G", "chr19:22614854:A:G")
index <- eur.df[eur.df$SNP %in% index,]

lowest_p_variants <- eur.df %>%
  group_by(CHR) %>%
  slice(which.min(P))

lowest_p_variants <- lowest_p_variants[!grepl("chr2\\:|chr8\\:|chr19\\:",lowest_p_variants$SNP),]
lowest_p_variants <- rbind.data.frame(lowest_p_variants, index)


df <- as.data.frame(eur.df$SNP)
colnames(df) <- "SNP"
# write.table(df, "EUR_5e-06_variants_clumping.txt", col.names = T, quote = F, row.names = F)


## Add EUR frequency
EUR.freq <- read.table("EUR.Alele1_freq_output_All.chr", header = T)
table(eur.df$SNP %in% EUR.freq$SNP )
eur.df$EUR.freq <- EUR.freq$MAF[match(eur.df$SNP, EUR.freq$SNP)]

## Add AFR frequency
AFR.freq <- read.table("EUR.Alele1_freq_output_All_in_replication_AFR.chr", header = T)
table(eur.df$SNP %in% AFR.freq$SNP )
eur.df$AFR.freq <- AFR.freq$MAF[match(eur.df$SNP, AFR.freq$SNP)]

## Add LD from Clumping
EUR.ld.clumping <- read.table("EUR_LD_clumped_variants_chr_ALL.ld", header = T)

eur.df$ld_with_index_var <- NULL
eur.df$index_var <- NULL
eur.df$index.var.pos <- NULL
for (chr in 1:22){
  print(chr)
  CHR= paste0("chr", chr, ":")
  index.var <- lowest_p_variants$SNP[grepl(CHR, lowest_p_variants$SNP)]
  index.var.pos <- lowest_p_variants$BP[grepl(CHR, lowest_p_variants$SNP)]
  eur.df$index_var[eur.df$CHR == chr] <- index.var
  eur.df$index.var.pos[eur.df$CHR == chr] <- index.var.pos
  EUR.ld.clumping.tmp <- EUR.ld.clumping[EUR.ld.clumping$CHR_A %in% chr,]
  if (nrow(EUR.ld.clumping.tmp) == 0) {
    next
  }
  
  # Swap SNP_A and SNP_B based on index.var to keep index always on the right column
  EUR.ld.clumping.tmp <- EUR.ld.clumping.tmp %>%
    mutate(
      temp_SNP = if_else(SNP_A == index.var, SNP_B, SNP_A),
      SNP_B = if_else(SNP_A == index.var, index.var, SNP_B),
      SNP_A = temp_SNP
    ) %>%
    select(-temp_SNP)
  EUR.ld.clumping.tmp <- EUR.ld.clumping.tmp[EUR.ld.clumping.tmp$SNP_B %in% index.var,]
  if (nrow(EUR.ld.clumping.tmp) == 0) {
    next
  }
  for (j in 1:nrow(EUR.ld.clumping.tmp)){
    other.var <- EUR.ld.clumping.tmp$SNP_A[j]
    eur.df$ld_with_index_var[eur.df$SNP == other.var] <- EUR.ld.clumping.tmp$R2[EUR.ld.clumping.tmp$SNP_A == other.var & EUR.ld.clumping.tmp$SNP_B == index.var]
  }
}

eur.df$distance_from_index_BP <- eur.df$BP - eur.df$index.var.pos

write.table(eur.df, "EUR_clumping_results.txt", col.names = T, quote = F, row.names = F)




## Clean GWAS tables
library(readxl)
library(data.table)
## read meta analysis file
metanalysis_data <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega/v2_with_preQC_VCF/T2D_Meta_analysis_SJLIFE_EUR_and_SJLIFE_AFR_fixed_1.tbl")

## Metasoft RE2
metasoft.res <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/METASOFT/metasoft_res_edited", sep = "\t", header=T)
header <- c("RSID", "STUDY", "PVALUE_FE", "BETA_FE", "STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2")
# header <- c("RSID", "STUDY", "PVALUE_FE", "BETA_FE", "STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE", "PVALUES_OF_STUDIES", "MVALUES_OF_STUDIES")

length(header)
metasoft.res <- metasoft.res[,1:11]
colnames(metasoft.res) <- header 
head(metasoft.res)


metanalysis_data.saved <- metanalysis_data
metanalysis_data.saved$Allele1 <- toupper(metanalysis_data.saved$Allele1)
metanalysis_data.saved$Allele2 <- toupper(metanalysis_data.saved$Allele2)

## FUnction to format P
format_P <- function(value) {
  if (value >= 0.1) {
    return(sprintf("%.2f", value))
  } else if (value >= 0.01) {
    return(sprintf("%.3f", value))
  } else {
    formatted <- formatC(value, format = "e", digits = 1)
    # # Replace "e" with "x10^"
    formatted <- (sub("e-0", "e-", formatted))
    return(sub("e", "x10^", formatted))
  }
}

##############
## AFR GWAS ##
##############
AFR.gwas <- read_excel("EDIT_TOP.AFR.only.with.P.5e-06.and.results.from.meta_cindy_cleaned.xlsx")


metanalysis_data <- metanalysis_data.saved[match(AFR.gwas$MarkerName, metanalysis_data.saved$MarkerName),]
metanalysis_data <- cbind.data.frame(metanalysis_data, AFR.gwas[match(metanalysis_data$MarkerName, AFR.gwas$MarkerName), c("A1", "A2")])
metanalysis_data$Effect[metanalysis_data$Allele1==metanalysis_data$A2] <- metanalysis_data$Effect[metanalysis_data$Allele1==metanalysis_data$A2]*-1

metanalysis_data$OddsRatio <- round(exp(metanalysis_data$Effect),2)

# Calculate the upper and lower confidence intervals for Odds Ratio
metanalysis_data$OR_Lower_CI <- round(exp(metanalysis_data$Effect - 1.96 * metanalysis_data$StdErr),2)
metanalysis_data$OR_Upper_CI <- round(exp(metanalysis_data$Effect + 1.96 * metanalysis_data$StdErr),2)

metanalysis_data$meta_OR <- paste0(metanalysis_data$OddsRatio, " (", metanalysis_data$OR_Lower_CI, "-", metanalysis_data$OR_Upper_CI, ")")


AFR.gwas <- cbind.data.frame(AFR.gwas, metanalysis_data[match(AFR.gwas$MarkerName, metanalysis_data$MarkerName),c("MarkerName", "meta_OR", "P-value", "HetPVal")])
AFR.gwas <- AFR.gwas[order(AFR.gwas$CHR), ]

AFR.gwas$NEW_OR_AFR <- paste0(round(AFR.gwas$OR_AFR,2), " (", round(AFR.gwas$L95_AFR,2), "-", round(AFR.gwas$U95_AFR,2), ")")
AFR.gwas$NEW_OR_EUR <- paste0(round(AFR.gwas$OR_EUR,2), " (", round(AFR.gwas$L95_EUR,2), "-", round(AFR.gwas$U95_EUR,2), ")")

AFR.gwas$ld_R2_with_index_var <- round(as.numeric(AFR.gwas$ld_R2_with_index_var),2)

## RE2
AFR.gwas <- cbind.data.frame(AFR.gwas, metasoft.res[match(AFR.gwas$MarkerName, metasoft.res$RSID), c("PVALUE_RE2", "STAT1_RE2",  "STAT2_RE2")])

## Reorder
AFR.gwas <- AFR.gwas[order(AFR.gwas$CHR, AFR.gwas$BP), ]

AFR.gwas$EAF_in_AFR <- round(AFR.gwas$EAF_in_AFR, 2)
AFR.gwas$EAF_in_EUR <- round(AFR.gwas$EAF_in_EUR, 2)

## Format P
AFR.gwas$P_AFR_cleaned <- sapply(AFR.gwas$P_AFR, format_P)
AFR.gwas$P_EUR_cleaned <- sapply(AFR.gwas$P_EUR, format_P)
AFR.gwas$HetPVal_cleaned <- sapply(AFR.gwas$HetPVal, format_P)
AFR.gwas$P_meta <- AFR.gwas$`P-value`
AFR.gwas$P_meta_cleaned <- sapply(AFR.gwas$P_meta, format_P)
AFR.gwas$P_meta_RE2 <- AFR.gwas$P_meta
AFR.gwas$P_meta_RE2[AFR.gwas$HetPVal<0.05] <- AFR.gwas$PVALUE_RE2[AFR.gwas$HetPVal<0.05]
AFR.gwas$P_meta_het_YN <- ifelse(AFR.gwas$HetPVal<0.05, "Yes", "No")
AFR.gwas$meta_OR_RE2_cleaned <- AFR.gwas$meta_OR
AFR.gwas$meta_OR_RE2_cleaned[AFR.gwas$HetPVal<0.05] <- ""
AFR.gwas$P_meta_RE2_cleaned <- sapply(AFR.gwas$P_meta_RE2, format_P)

write.table(AFR.gwas, "AFR_gwas_cleaned_results_v2.txt", col.names = T, quote = F, row.names = F, sep = "\t")

##############
## EUR GWAS ##
##############
EUR.gwas <- read_excel("EDIT_chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA.P5e-06_Yadav_Cindy_cleaned.xlsx")
EUR.gwas$MarkerName <- EUR.gwas$SNP

## Fix Meta beta
metanalysis_data <- metanalysis_data.saved[match(EUR.gwas$MarkerName, metanalysis_data.saved$MarkerName),]
metanalysis_data <- cbind.data.frame(metanalysis_data, EUR.gwas[match(metanalysis_data$MarkerName, EUR.gwas$MarkerName), c("A1", "A2")])
metanalysis_data$Effect[metanalysis_data$Allele1==metanalysis_data$A2] <- metanalysis_data$Effect[metanalysis_data$Allele1==metanalysis_data$A2]*-1

metanalysis_data$OddsRatio <- round(exp(metanalysis_data$Effect),2)
# Calculate the upper and lower confidence intervals for Odds Ratio
metanalysis_data$OR_Lower_CI <- round(exp(metanalysis_data$Effect - 1.96 * metanalysis_data$StdErr),2)
metanalysis_data$OR_Upper_CI <- round(exp(metanalysis_data$Effect + 1.96 * metanalysis_data$StdErr),2)

metanalysis_data$meta_OR <- paste0(metanalysis_data$OddsRatio, " (", metanalysis_data$OR_Lower_CI, "-", metanalysis_data$OR_Upper_CI, ")")

EUR.gwas <- cbind.data.frame(EUR.gwas, metanalysis_data[match(EUR.gwas$SNP, metanalysis_data$MarkerName),c("MarkerName", "meta_OR", "P-value", "HetPVal")])


EUR.gwas$NEW_OR_EUR <- paste0(round(EUR.gwas$OR_EUR,2), " (", round(EUR.gwas$L95_EUR,2), "-", round(EUR.gwas$U95_EUR,2), ")")
EUR.gwas$NEW_OR_AFR <- paste0(round(EUR.gwas$OR_AFR,2), " (", round(EUR.gwas$L95_AFR,2), "-", round(EUR.gwas$U95_AFR,2), ")")

EUR.gwas$ld_R2_with_index_var <- round(as.numeric(EUR.gwas$ld_R2_with_index_var),2)


## RE2
EUR.gwas <- cbind.data.frame(EUR.gwas, metasoft.res[match(EUR.gwas$MarkerName, metasoft.res$RSID), c("PVALUE_RE2", "STAT1_RE2",  "STAT2_RE2")])

## Reorder
EUR.gwas <- EUR.gwas[order(EUR.gwas$CHR, EUR.gwas$BP), ]

EUR.gwas$EAF_in_EUR <- round(EUR.gwas$EAF_in_EUR, 2)
EUR.gwas$EAF_in_AFR <- round(EUR.gwas$EAF_in_AFR, 2)


## Format P
EUR.gwas$P_EUR_cleaned <- sapply(EUR.gwas$P_EUR, format_P)
EUR.gwas$P_AFR_cleaned <- sapply(EUR.gwas$P_AFR, format_P)
EUR.gwas$HetPVal_cleaned <- sapply(EUR.gwas$HetPVal, format_P)
EUR.gwas$P_meta <- EUR.gwas$`P-value`
EUR.gwas$P_meta_cleaned <- sapply(EUR.gwas$P_meta, format_P)
EUR.gwas$P_meta_RE2 <- EUR.gwas$P_meta
EUR.gwas$P_meta_RE2[EUR.gwas$HetPVal<0.05] <- EUR.gwas$PVALUE_RE2[EUR.gwas$HetPVal<0.05]
EUR.gwas$P_meta_het_YN <- ifelse(EUR.gwas$HetPVal<0.05, "Yes", "No")
EUR.gwas$meta_OR_RE2_cleaned <- EUR.gwas$meta_OR
EUR.gwas$meta_OR_RE2_cleaned[EUR.gwas$HetPVal<0.05] <- ""
EUR.gwas$P_meta_RE2_cleaned <- sapply(EUR.gwas$P_meta_RE2, format_P)

write.table(EUR.gwas, "EUR_gwas_cleaned_results_v2.txt", col.names = T, quote = F, row.names = F, sep = "\t")

rm(metanalysis_data)
