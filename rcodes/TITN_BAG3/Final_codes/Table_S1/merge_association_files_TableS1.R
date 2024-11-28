combined.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples.bim", header = F, sep = "")
combined.bim$KEY1 <- paste0(combined.bim$V1, ":", combined.bim$V4)
sum(duplicated(combined.bim$KEY1))
## Harmonize
combined.bim$KEY <- paste0("chr", combined.bim$V1, ":", combined.bim$V4, ":", combined.bim$V6, ":", combined.bim$V5)
sum(duplicated(combined.bim$KEY))

keys <- combined.bim$KEY[duplicated(combined.bim$KEY)]
harmonize <- combined.bim[combined.bim$KEY %in% keys,]
harmonize <- harmonize[c("V2", "KEY")]
# write.table(harmonize, "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/harmonize.txt", sep = "\t", col.names = F, row.names = F, quote = F)

##############
## Combined ##
##############
combined <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples_results_kendrick.assoc.logistic", header = T, sep = "")
dim(combined)
combined.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples_updated.bim", header = F, sep = "")
combined.bim$KEY1 <- paste0(combined.bim$V1, ":", combined.bim$V4)
sum(duplicated(combined.bim$KEY1))
combined$REF <- combined.bim$V6[match(combined$SNP, combined.bim$V2)]
combined$ALT <- combined.bim$V5[match(combined$SNP, combined.bim$V2)]
table(combined$ALT == combined$A1)
swapped <- combined$ALT != combined$A1
# Swap ALT and REF for these rows
combined[swapped, c("ALT", "REF")] <- combined[swapped, c("REF", "ALT")]
# Adjust OR, L95, and U95 for swapped alleles
combined[swapped, c("OR", "L95", "U95")] <- 1 / combined[swapped, c("OR", "U95", "L95")]
table(combined$ALT == combined$A1)
combined$OR_95CI <- paste0(round(combined$OR,2), " (", round(combined$L95, 2), "-", round(combined$U95,2), ")")
combined$KEY <- paste0("chr", combined$CHR, ":", combined$BP, ":", combined$REF, ":", combined$A1)  
combined <- combined[c("SNP", "CHR", "BP", "A1", "REF", "OR_95CI", "P", "KEY")]


############
## SJLIFE ##
############
sjlife <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_results_kendrick.assoc.logistic", header = T, sep = "")
dim(sjlife)
sjlife.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_to_concat_updated.bim", header = F, sep = "")
sjlife.bim$KEY1 <- paste0(sjlife.bim$V1, ":", sjlife.bim$V4)
sum(duplicated(sjlife.bim$KEY1))

sjlife$REF <- sjlife.bim$V6[match(sjlife$SNP, sjlife.bim$V2)]
sjlife$ALT <- sjlife.bim$V5[match(sjlife$SNP, sjlife.bim$V2)]
table(sjlife$ALT == sjlife$A1)
swapped <- sjlife$ALT != sjlife$A1
# Swap ALT and REF for these rows
sjlife[swapped, c("ALT", "REF")] <- sjlife[swapped, c("REF", "ALT")]
# Adjust OR, L95, and U95 for swapped alleles
sjlife[swapped, c("OR", "L95", "U95")] <- 1 / sjlife[swapped, c("OR", "U95", "L95")]
sjlife$OR_95CI <- paste0(round(sjlife$OR,2), " (", round(sjlife$L95, 2), "-", round(sjlife$U95,2), ")")
sjlife$KEY <- paste0("chr", sjlife$CHR, ":", sjlife$BP, ":", sjlife$REF, ":", sjlife$A1)  
sjlife <- sjlife[c("OR_95CI", "P", "KEY", "A1")]

# ##############
# ## CCSS_org ##
# ##############
# ccss_org <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_org.assoc.logistic", header = T, sep = "")
# dim(ccss_org)
# ccss_org.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_org_to_concat_updated.bim", header = F, sep = "")
# ccss_org.bim$KEY1 <- paste0(ccss_org.bim$V1, ":", ccss_org.bim$V4)
# sum(duplicated(ccss_org.bim$KEY1))
# ccss_org$REF <- ccss_org.bim$V6[match(ccss_org$SNP, ccss_org.bim$V2)]
# ccss_org$ALT <- ccss_org.bim$V5[match(ccss_org$SNP, ccss_org.bim$V2)]
# table(ccss_org$ALT == ccss_org$A1)
# ccss_org[ccss_org$ALT != ccss_org$A1, c("ALT", "REF")] <- ccss_org[ccss_org$ALT != ccss_org$A1, c("REF", "ALT")]
# ccss_org$OR_95CI <- paste0(round(ccss_org$OR,2), " (", round(ccss_org$L95, 2), "-", round(ccss_org$U95,2), ")")
# ccss_org$KEY <- paste0("chr", ccss_org$CHR, ":", ccss_org$BP, ":", ccss_org$REF, ":", ccss_org$A1)
# ccss_org <- ccss_org[c("OR_95CI", "P", "KEY")]
# 
# ##############
# ## CCSS_exp ##
# ##############
# ccss_exp <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_exp.assoc.logistic", header = T, sep = "")
# dim(ccss_exp)
# ccss_exp.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_exp_to_concat_updated.bim", header = F, sep = "")
# ccss_exp.bim$KEY1 <- paste0(ccss_exp.bim$V1, ":", ccss_exp.bim$V4)
# sum(duplicated(ccss_exp.bim$KEY1))
# ccss_exp$REF <- ccss_exp.bim$V6[match(ccss_exp$SNP, ccss_exp.bim$V2)]
# ccss_exp$ALT <- ccss_exp.bim$V5[match(ccss_exp$SNP, ccss_exp.bim$V2)]
# table(ccss_exp$ALT == ccss_exp$A1)
# ccss_exp[ccss_exp$ALT != ccss_exp$A1, c("ALT", "REF")] <- ccss_exp[ccss_exp$ALT != ccss_exp$A1, c("REF", "ALT")]
# ccss_exp$OR_95CI <- paste0(round(ccss_exp$OR,2), " (", round(ccss_exp$L95, 2), "-", round(ccss_exp$U95,2), ")")
# ccss_exp$KEY <- paste0("chr", ccss_exp$CHR, ":", ccss_exp$BP, ":", ccss_exp$REF, ":", ccss_exp$A1)
# ccss_exp <- ccss_exp[c("OR_95CI", "P", "KEY")]

###################
## CCSS_combined ##
###################
ccss_org_ccss_exp <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_org_ccss_exp_kendrick.assoc.logistic", header = T, sep = "")
dim(ccss_org_ccss_exp)
ccss_org_ccss_exp.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/merged_ccss.bim", header = F, sep = "")
ccss_org_ccss_exp.bim$KEY1 <- paste0(ccss_org_ccss_exp.bim$V1, ":", ccss_org_ccss_exp.bim$V4)
sum(duplicated(ccss_org_ccss_exp.bim$KEY1))
ccss_org_ccss_exp$REF <- ccss_org_ccss_exp.bim$V6[match(ccss_org_ccss_exp$SNP, ccss_org_ccss_exp.bim$V2)]
ccss_org_ccss_exp$ALT <- ccss_org_ccss_exp.bim$V5[match(ccss_org_ccss_exp$SNP, ccss_org_ccss_exp.bim$V2)]
table(ccss_org_ccss_exp$ALT != ccss_org_ccss_exp$A1)
swapped <- ccss_org_ccss_exp$ALT != ccss_org_ccss_exp$A1
# Swap ALT and REF for these rows
ccss_org_ccss_exp[swapped, c("ALT", "REF")] <- ccss_org_ccss_exp[swapped, c("REF", "ALT")]
# Adjust OR, L95, and U95 for swapped alleles
ccss_org_ccss_exp[swapped, c("OR", "L95", "U95")] <- 1 / ccss_org_ccss_exp[swapped, c("OR", "U95", "L95")]

ccss_org_ccss_exp$OR_95CI <- paste0(round(ccss_org_ccss_exp$OR,2), " (", round(ccss_org_ccss_exp$L95, 2), "-", round(ccss_org_ccss_exp$U95,2), ")")
ccss_org_ccss_exp$KEY <- paste0("chr", ccss_org_ccss_exp$CHR, ":", ccss_org_ccss_exp$BP, ":", ccss_org_ccss_exp$REF, ":", ccss_org_ccss_exp$A1)  
ccss_org_ccss_exp <- ccss_org_ccss_exp[c("OR_95CI", "P", "KEY", "A1")]



library(purrr)

# dfs <- list(combined, sjlife, ccss_org, ccss_exp, ccss_org_ccss_exp, sjlife_afr)
dfs <- list(combined, sjlife, ccss_org_ccss_exp)
final_merged_df <- reduce(dfs, merge, by = "KEY")
# colnames(final_merged_df) <- c("SNP", "KEY", "CHR", "BP", "A1", "REF", "OR_95CI_combined", "P_combined", 
#                                "OR_95CI_sjlife", "P_sjlife", 
#                                "OR_95CI_ccss_org", "P_ccss_org", 
#                                "OR_95CI_ccss_exp", "P_ccss_exp", 
#                                "OR_95CI_ccss_org_exp", "P_ccss_org_exp")

colnames(final_merged_df) <- c("SNP", "KEY", "CHR", "BP", "A1", "REF", "OR_95CI_combined", "P_combined", 
                               "OR_95CI_sjlife", "P_sjlife", "A1_sjlife", 
                               "OR_95CI_ccss_org_exp", "P_ccss_org_exp", "A1_ccss_org_exp")


freq.df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples_updated_freq_out_kendrick.frq", header = T)
dim(freq.df)
final_merged_df$combinedMAF <- freq.df$MAF[match(final_merged_df$SNP, freq.df$SNP)]
# final_merged_df$freqA1 <- freq.df$A1[match(final_merged_df$SNP, freq.df$SNP)]

## Filter for MAF > 0.05
final_merged_df <- final_merged_df[final_merged_df$combinedMAF > 0.05,]

###########################
## read annoatation file ##
###########################
library(data.table)

TTN <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/chr2_178_list_ALL_08_27_2024.txt", header = T)
TTN <- cbind.data.frame(CHR=TTN$CHROM, POS=TTN$POS, SNP=TTN$ID, REF=TTN$REF, ALT=TTN$ALT, 
                        GENE=TTN$`ANN[*].GENE`, EFFECT=TTN$`ANN[*].EFFECT`, AF=TTN$AF, AF_nfe=TTN$AF_nfe, AF_afr=TTN$AF_afr)

BAG3 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/chr10_119_list_ALL_08_27_2024.txt", header = T)
BAG3 <- cbind.data.frame(CHR=BAG3$CHROM, POS=BAG3$POS, SNP=BAG3$ID, REF=BAG3$REF, ALT=BAG3$ALT, 
                         GENE=BAG3$`ANN[*].GENE`, EFFECT=BAG3$`ANN[*].EFFECT`, AF=BAG3$AF, AF_nfe=BAG3$AF_nfe, AF_afr=BAG3$AF_afr)

annotation <- rbind.data.frame(TTN, BAG3)
# annotation$SNPID <- sub(";rs[0-9]+", "", annotation$SNP)
# annotation$RSID <- sub(".*;rs([0-9]+)", "rs\\1", annotation$SNP)
annotation$SNPID <- sub(";.*", "", annotation$SNP)
annotation$RSID <- sub(".*?;([^;]*rs[^;]*).*", "\\1", annotation$SNP)
annotation$RSID[!grepl("rs", annotation$rsID)] <- ""

## extract missense
annotation <- annotation[grepl("missense", annotation$EFFECT, ignore.case = T),]
# sum(final_merged_df$KEY %in% annotation$SNPID)

## Match by both REF:ALT and ALT:REF
annotation$KEY1 <- paste0(annotation$CHR, ":", annotation$POS, ":", annotation$REF, ":", annotation$ALT)
annotation$KEY2 <- paste0(annotation$CHR, ":", annotation$POS, ":", annotation$ALT, ":",annotation$REF)

annotation$annotationKEY <- final_merged_df$KEY[match(annotation$KEY1, final_merged_df$KEY)]
# If there are NAs, try to match them with KEY2
missing_indices <- is.na(annotation$annotationKEY)
annotation$annotationKEY[missing_indices] <- final_merged_df$KEY[match(annotation$KEY2[missing_indices], final_merged_df$KEY)]
annotation <- annotation[!is.na(annotation$annotationKEY),]
cc <- final_merged_df[final_merged_df$KEY %in% annotation$annotationKEY,]
dim(cc)
# 40 15

cc$variation <- annotation$EFFECT[match(cc$KEY, annotation$annotationKEY)]
cc$Gene <- annotation$GENE[match(cc$KEY, annotation$annotationKEY)]
cc$rsID <- annotation$RSID[match(cc$KEY, annotation$annotationKEY)]
cc <- cc[order(cc$P_combined), ]
cc <- cc[grepl("TTN|BAG3", cc$Gene),]

################
## SJLIFE afr ##
################
sjlife_afr <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_results_afr_AN_kendrick.assoc.logistic", header = T, sep = "")
dim(sjlife_afr)
sjlife_afr.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_afr_to_concat_updated.bim", header = F, sep = "")
sjlife_afr.bim$KEY1 <- paste0(sjlife_afr.bim$V1, ":", sjlife_afr.bim$V4)
sum(duplicated(sjlife_afr.bim$KEY1))

afr.frq <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_results_afr_freq_AN_kendrick.frq", header = T)
table(sjlife_afr$SNP %in% afr.frq$SNP)
sjlife_afr$afr_MAF <- afr.frq$MAF[match(sjlife_afr$SNP, afr.frq$SNP)]
sjlife_afr$afr_freqA1 <- afr.frq$A1[match(sjlife_afr$SNP, afr.frq$SNP)]

sjlife_afr$REF <- sjlife_afr.bim$V6[match(sjlife_afr$SNP, sjlife_afr.bim$V2)]
sjlife_afr$ALT <- sjlife_afr.bim$V5[match(sjlife_afr$SNP, sjlife_afr.bim$V2)]
table(sjlife_afr$ALT != sjlife_afr$A1)
# swapped <- sjlife_afr$ALT != sjlife_afr$A1
# sjlife_afr[swapped,]
# # Swap ALT and REF for these rows
# sjlife_afr[swapped, c("ALT", "REF")] <- sjlife_afr[swapped, c("REF", "ALT")]
# # Adjust OR, L95, and U95 for swapped alleles
# sjlife_afr[swapped, c("OR", "L95", "U95")] <- 1 / sjlife_afr[swapped, c("OR", "U95", "L95")]
# sjlife_afr$OR_95CI_afr <- paste0(round(sjlife_afr$OR,2), " (", round(sjlife_afr$L95, 2), "-", round(sjlife_afr$U95,2), ")")
sjlife_afr$KEY <- paste0("chr", sjlife_afr$CHR, ":", sjlife_afr$BP, ":", sjlife_afr$REF, ":", sjlife_afr$A1)  
# sjlife_afr <- sjlife_afr[c("OR_95CI_afr", "P", "KEY", "afr_MAF", "SNP", "A1")]
# cc$OR_95CI_afr <- sjlife_afr$OR_95CI_afr[match(cc$KEY, sjlife_afr$SNP)]
cc$OR_afr <- sjlife_afr$OR[match(cc$KEY, sjlife_afr$SNP)]
cc$L95_afr <- sjlife_afr$L95[match(cc$KEY, sjlife_afr$SNP)]
cc$U95_afr <- sjlife_afr$U95[match(cc$KEY, sjlife_afr$SNP)]
cc$P_afr <- sjlife_afr$P[match(cc$KEY, sjlife_afr$SNP)]
cc$afr_MAF <- sjlife_afr$afr_MAF[match(cc$KEY, sjlife_afr$SNP)]
cc$afr_freqA1 <- sjlife_afr$afr_freqA1[match(cc$KEY, sjlife_afr$SNP)]
cc$A1_sjlife_afr <- sjlife_afr$A1[match(cc$KEY, sjlife_afr$SNP)]
# cc$ALT_sjlife_afr <- sjlife_afr$ALT[match(cc$KEY, sjlife_afr$SNP)]

table(cc$A1 != cc$A1_sjlife_afr)
swapped <- cc$A1 != cc$A1_sjlife_afr
# Swap ALT and REF for these rows
# cc[swapped, c("ALT", "REF")] <- cc[swapped, c("REF", "ALT")]
# Adjust OR, L95, and U95 for swapped alleles
cc[swapped, c("OR_afr", "L95_afr", "U95_afr")] <- 1 / cc[swapped, c("OR_afr", "U95_afr", "L95_afr")]
cc$OR_95CI_afr <- paste0(round(cc$OR_afr,2), " (", round(cc$L95_afr, 2), "-", round(cc$U95_afr,2), ")")


#############################################
## Create meta analysis summary statistics ##
#############################################
# We need MarkerName Effect StdErr P N
## SJLIFE
sjlife <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_results_kendrick.assoc.logistic", header = T, sep = "")
sjlife$BETA <- log(sjlife$OR)
SJLIFE.metal <- cbind.data.frame(SNP=cc$SNP, cc$OR_95CI_sjlife, P=cc$P_sjlife)
SJLIFE.metal$N <- sjlife$NMISS[match(SJLIFE.metal$SNP, sjlife$SNP)]
SJLIFE.metal$SE <- sjlife$SE[match(SJLIFE.metal$SNP, sjlife$SNP)]
SJLIFE.metal$BETA <- sjlife$BETA[match(SJLIFE.metal$SNP, sjlife$SNP)]
SJLIFE.metal <- cbind.data.frame(MarkerName=SJLIFE.metal$SNP, Effect=SJLIFE.metal$BETA, StdErr=SJLIFE.metal$SE, P=SJLIFE.metal$P, N=SJLIFE.metal$N)
write.table(SJLIFE.metal, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/common_variant_analysis_SJLIFE.txt", col.names  = T, sep = " ", row.names = F, quote = F)

## CCSS
CCSS <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_org_ccss_exp_kendrick.assoc.logistic", header = T, sep = "")
CCSS$BETA <- log(CCSS$OR)
CCSS.metal <- cbind.data.frame(SNP=cc$SNP, cc$OR_95CI_ccss_org_exp, P=cc$P_ccss_org_exp)
CCSS.metal$N <- CCSS$NMISS[match(CCSS.metal$SNP, CCSS$SNP)]
CCSS.metal$SE <- CCSS$SE[match(CCSS.metal$SNP, CCSS$SNP)]
CCSS.metal$BETA <- CCSS$BETA[match(CCSS.metal$SNP, CCSS$SNP)]
CCSS.metal <- cbind.data.frame(MarkerName=CCSS.metal$SNP, Effect=CCSS.metal$BETA, StdErr=CCSS.metal$SE, P=CCSS.metal$P, N=CCSS.metal$N)
write.table(CCSS.metal, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/common_variant_analysis_CCSS.txt", col.names  = T, sep = " ", row.names = F, quote = F)

## SJLIFE AFR
sjlife.afr <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_results_afr_AN_kendrick.assoc.logistic", header = T, sep = "")
sjlife.afr$BETA <- log(sjlife.afr$OR)
SJLIFE.AFR.metal <- cbind.data.frame(SNP=cc$SNP, cc$OR_95CI_afr, P=cc$P_afr)
SJLIFE.AFR.metal$N <- sjlife.afr$NMISS[match(SJLIFE.AFR.metal$SNP, sjlife.afr$SNP)]
SJLIFE.AFR.metal$SE <- sjlife.afr$SE[match(SJLIFE.AFR.metal$SNP, sjlife.afr$SNP)]
SJLIFE.AFR.metal$BETA <- sjlife.afr$BETA[match(SJLIFE.AFR.metal$SNP, sjlife.afr$SNP)]
SJLIFE.AFR.metal <- cbind.data.frame(MarkerName=SJLIFE.AFR.metal$SNP, Effect=SJLIFE.AFR.metal$BETA, StdErr=SJLIFE.AFR.metal$SE, P=SJLIFE.AFR.metal$P, N=SJLIFE.AFR.metal$N)
write.table(SJLIFE.AFR.metal, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/common_variant_analysis_SJLIFE.AFR.txt", col.names  = T, sep = " ", row.names = F, quote = F)



## After running metal metal_scripts_common_variants.txt in /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis:
data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/common_variant_Meta_analysis_sjlife_ccss_fixed_stratified_analysis_1.tbl", header = TRUE, sep = "\t")

# Calculate OR_CI as effect size and standard error
data$OR_CI <- paste(round(exp(data$Effect), 3), " (", round(exp(data$Effect - 1.96 * data$StdErr), 3), "-", round(exp(data$Effect + 1.96 * data$StdErr), 3), ")", sep = "")

# Select and display relevant columns
result <- data[, c("MarkerName", "P.value", "OR_CI")]
# sor the results in the same order
cc$metaOR <- result$OR_CI[match(cc$SNP, result$MarkerName)]
cc$metaP <- result$P.value[match(cc$SNP, result$MarkerName)]

# sort based on existing table
snp_vector <- c("rs744426", "rs3731746", "rs2303838", "rs2042996", "rs9808377", 
                "rs3829747", "rs1001238", "rs16866406", "rs2288569", "rs2042995", 
                "rs2627043", "rs3731749", "rs6723526", "rs16866465", "rs12693166", 
                "rs12693164", "rs13390491", "rs16866538", "rs35813871", "rs10497520", 
                "rs12463674", "rs36051007", "rs2244492", "rs72648998", "rs2163008", 
                "rs72648907", "rs2291310", "rs2291311", "rs7585334", "rs2627037", 
                "rs922984", "rs3858340", "rs3829746", "rs2234962")


cc_sorted <- cc[match(snp_vector, cc$rsID), ]

cc_sorted <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/common_missense_variants_from_kendricks_pheno.txt", header = T, sep = "\t")

## Add frequency for both SJLIFE and CCSS separately
freq.df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_to_concat_updated_freq_out_kendrick.frq", header = T)
cc_sorted$SJLIFE.EUR.MAF <- freq.df$MAF[match(cc_sorted$SNP, freq.df$SNP)]  

freq.df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/merged_ccss_freq_out_kendrick.frq", header = T)
cc_sorted$CCSS.EUR.MAF <- freq.df$MAF[match(cc_sorted$SNP, freq.df$SNP)]  

## Add metal analysis results with SJLIFE, CCSS and AFR
## After running metal metal_scripts_common_variants_ccss_sjlife_afr.txt in /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis:
data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/common_variant_Meta_analysis_sjlife_ccss_afr_fixed_stratified_analysis_1.tbl", header = TRUE, sep = "\t")

# Calculate OR_CI as effect size and standard error
data$OR_CI <- paste(round(exp(data$Effect), 3), " (", round(exp(data$Effect - 1.96 * data$StdErr), 3), "-", round(exp(data$Effect + 1.96 * data$StdErr), 3), ")", sep = "")

# Select and display relevant columns
result <- data[, c("MarkerName", "P.value", "OR_CI")]
# sor the results in the same order
cc_sorted$metaOR.withAFR <- result$OR_CI[match(cc_sorted$SNP, result$MarkerName)]
cc_sorted$metaP.with.AFR <- result$P.value[match(cc_sorted$SNP, result$MarkerName)]



write.table(cc_sorted, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/common_missense_variants_from_kendricks_pheno_with_meta_analysis_with_SJLIFE_CCSS_and_AFR.txt", col.names = T, sep = "\t", quote = F, row.names = F)

# reference <- read.table(text="chr	pos	EA	NEA
# 2	178633315	C	T
# 2	178580212	T	C
# 2	178566270	A	G
# 2	178571293	A	G
# 2	178586693	A	G
# 2	178556967	G	A
# 2	178599800	C	T
# 2	178532834	T	C
# 2	178541464	T	C
# 2	178592420	A	G
# 2	178562809	C	T
# 2	178693639	C	T
# 2	178593864	T	C
# 2	178718769	G	T
# 2	178722403	G	C
# 2	178714366	C	T
# 2	178717600	T	C
# 2	178785681	A	G
# 2	178795185	A	G
# 2	178717810	T	G
# 2	178567458	G	A
# 2	178681132	T	C
# 2	178780128	T	C
# 2	178689578	T	C
# 2	178759031	C	T
# 2	178764734	C	T
# 2	178756224	C	T
# 2	178751160	T	C
# 2	178741811	A	G
# 2	178663651	T	C
# 2	178747656	T	C
# 2	178710784	T	C
# 10	119670121	C	T
# 10	119676774	T	C", header = T)
# 
# reference$KEY <- paste0("chr", reference$chr,":", reference$pos)
