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

## read meta analysis file
metanalysis_data <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega/v2_with_preQC_VCF/T2D_Meta_analysis_SJLIFE_EUR_and_SJLIFE_AFR_fixed_1.tbl")

metanalysis_data$OddsRatio <- round(exp(metanalysis_data$Effect),2)

# Calculate the upper and lower confidence intervals for Odds Ratio
metanalysis_data$OR_Lower_CI <- round(exp(metanalysis_data$Effect - 1.96 * metanalysis_data$StdErr),2)
metanalysis_data$OR_Upper_CI <- round(exp(metanalysis_data$Effect + 1.96 * metanalysis_data$StdErr),2)

metanalysis_data$meta_OR <- paste0(metanalysis_data$OddsRatio, " (", metanalysis_data$OR_Lower_CI, "-", metanalysis_data$OR_Upper_CI, ")")

## AFR GWAS
AFR.gwas <- read_excel("EDIT_TOP.AFR.only.with.P.5e-06.and.results.from.meta_cindy_cleaned.xlsx")
AFR.gwas <- cbind.data.frame(AFR.gwas, metanalysis_data[match(AFR.gwas$MarkerName, metanalysis_data$MarkerName),c("MarkerName", "meta_OR", "P-value", "HetPVal")])
AFR.gwas <- AFR.gwas[order(AFR.gwas$CHR), ]

AFR.gwas$NEW_OR_AFR <- paste0(AFR.gwas$OR_AFR, " (", AFR.gwas$L95_AFR, "-", AFR.gwas$U95_AFR, ")")
AFR.gwas$NEW_OR_EUR <- paste0(AFR.gwas$OR_EUR, " (", AFR.gwas$L95_EUR, "-", AFR.gwas$U95_EUR, ")")

AFR.gwas$ld_R2_with_index_var <- as.numeric(AFR.gwas$ld_R2_with_index_var)


write.table(AFR.gwas, "AFR_gwas_cleaned_results.txt", col.names = T, quote = F, row.names = F, sep = "\t")

## EUR GWAS
EUR.gwas <- read_excel("EDIT_chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA.P5e-06_Yadav_Cindy_cleaned.xlsx")
EUR.gwas <- cbind.data.frame(EUR.gwas, metanalysis_data[match(EUR.gwas$SNP, metanalysis_data$MarkerName),c("MarkerName", "meta_OR", "P-value", "HetPVal")])
EUR.gwas <- EUR.gwas[order(EUR.gwas$CHR), ]

EUR.gwas$NEW_OR_EUR <- paste0(EUR.gwas$OR_EUR, " (", EUR.gwas$L95_EUR, "-", EUR.gwas$U95_EUR, ")")
EUR.gwas$NEW_OR_AFR <- paste0(EUR.gwas$OR_AFR, " (", EUR.gwas$L95_AFR, "-", EUR.gwas$U95_AFR, ")")

EUR.gwas$ld_R2_with_index_var <- as.numeric(EUR.gwas$ld_R2_with_index_var)




write.table(EUR.gwas, "EUR_gwas_cleaned_results.txt", col.names = T, quote = F, row.names = F, sep = "\t")

rm(metanalysis_data)
