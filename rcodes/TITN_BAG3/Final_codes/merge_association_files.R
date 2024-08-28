combined.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples.bim", header = F, sep = "")
combined.bim$KEY1 <- paste0(combined.bim$V1, ":", combined.bim$V4)
sum(duplicated(combined.bim$KEY1))
## Harmonize
combined.bim$KEY <- paste0("chr", combined.bim$V1, ":", combined.bim$V4, ":", combined.bim$V6, ":", combined.bim$V5)
sum(duplicated(combined.bim$KEY))

keys <- combined.bim$KEY[duplicated(combined.bim$KEY)]
harmonize <- combined.bim[combined.bim$KEY %in% keys,]
harmonize <- harmonize[c("V2", "KEY")]
write.table(harmonize, "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/harmonize.txt", sep = "\t", col.names = F, row.names = F, quote = F)


# Combined
combined <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples_results.assoc.logistic", header = T, sep = "")
dim(combined)
combined.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples_updated.bim", header = F, sep = "")
combined.bim$KEY1 <- paste0(combined.bim$V1, ":", combined.bim$V4)
sum(duplicated(combined.bim$KEY1))
combined$REF <- combined.bim$V6[match(combined$SNP, combined.bim$V2)]
combined$ALT <- combined.bim$V5[match(combined$SNP, combined.bim$V2)]
table(combined$ALT == combined$A1)
combined[combined$ALT != combined$A1, c("ALT", "REF")] <- combined[combined$ALT != combined$A1, c("REF", "ALT")]
table(combined$ALT == combined$A1)
combined$OR_95CI <- paste0(round(combined$OR,2), " (", round(combined$L95, 2), "-", round(combined$U95,2), ")")
combined$KEY <- paste0("chr", combined$CHR, ":", combined$BP, ":", combined$REF, ":", combined$A1)  
combined <- combined[c("SNP", "CHR", "BP", "A1", "REF", "OR_95CI", "P", "KEY")]

# SJLIFE
sjlife <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_results.assoc.logistic", header = T, sep = "")
dim(sjlife)
sjlife.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/sjlife_to_concat_updated.bim", header = F, sep = "")
sjlife.bim$KEY1 <- paste0(sjlife.bim$V1, ":", sjlife.bim$V4)
sum(duplicated(sjlife.bim$KEY1))
sjlife$REF <- sjlife.bim$V6[match(sjlife$SNP, sjlife.bim$V2)]
sjlife$ALT <- sjlife.bim$V5[match(sjlife$SNP, sjlife.bim$V2)]
table(sjlife$ALT == sjlife$A1)
sjlife[sjlife$ALT != sjlife$A1, c("ALT", "REF")] <- sjlife[sjlife$ALT != sjlife$A1, c("REF", "ALT")]
sjlife$OR_95CI <- paste0(round(sjlife$OR,2), " (", round(sjlife$L95, 2), "-", round(sjlife$U95,2), ")")
sjlife$KEY <- paste0("chr", sjlife$CHR, ":", sjlife$BP, ":", sjlife$REF, ":", sjlife$A1)  
sjlife <- sjlife[c("OR_95CI", "P", "KEY")]

# CCSS_org
ccss_org <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_org.assoc.logistic", header = T, sep = "")
dim(ccss_org)
ccss_org.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_org_to_concat_updated.bim", header = F, sep = "")
ccss_org.bim$KEY1 <- paste0(ccss_org.bim$V1, ":", ccss_org.bim$V4)
sum(duplicated(ccss_org.bim$KEY1))
ccss_org$REF <- ccss_org.bim$V6[match(ccss_org$SNP, ccss_org.bim$V2)]
ccss_org$ALT <- ccss_org.bim$V5[match(ccss_org$SNP, ccss_org.bim$V2)]
table(ccss_org$ALT == ccss_org$A1)
ccss_org[ccss_org$ALT != ccss_org$A1, c("ALT", "REF")] <- ccss_org[ccss_org$ALT != ccss_org$A1, c("REF", "ALT")]
ccss_org$OR_95CI <- paste0(round(ccss_org$OR,2), " (", round(ccss_org$L95, 2), "-", round(ccss_org$U95,2), ")")
ccss_org$KEY <- paste0("chr", ccss_org$CHR, ":", ccss_org$BP, ":", ccss_org$REF, ":", ccss_org$A1)  
ccss_org <- ccss_org[c("OR_95CI", "P", "KEY")]

## CCSS_exp
ccss_exp <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_exp.assoc.logistic", header = T, sep = "")
dim(ccss_exp)
ccss_exp.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_exp_to_concat_updated.bim", header = F, sep = "")
ccss_exp.bim$KEY1 <- paste0(ccss_exp.bim$V1, ":", ccss_exp.bim$V4)
sum(duplicated(ccss_exp.bim$KEY1))
ccss_exp$REF <- ccss_exp.bim$V6[match(ccss_exp$SNP, ccss_exp.bim$V2)]
ccss_exp$ALT <- ccss_exp.bim$V5[match(ccss_exp$SNP, ccss_exp.bim$V2)]
table(ccss_exp$ALT == ccss_exp$A1)
ccss_exp[ccss_exp$ALT != ccss_exp$A1, c("ALT", "REF")] <- ccss_exp[ccss_exp$ALT != ccss_exp$A1, c("REF", "ALT")]
ccss_exp$OR_95CI <- paste0(round(ccss_exp$OR,2), " (", round(ccss_exp$L95, 2), "-", round(ccss_exp$U95,2), ")")
ccss_exp$KEY <- paste0("chr", ccss_exp$CHR, ":", ccss_exp$BP, ":", ccss_exp$REF, ":", ccss_exp$A1)  
ccss_exp <- ccss_exp[c("OR_95CI", "P", "KEY")]

## CCSS_combined
ccss_org_ccss_exp <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ccss_org_ccss_exp.assoc.logistic", header = T, sep = "")
dim(ccss_org_ccss_exp)
ccss_org_ccss_exp.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/merged_ccss.bim", header = F, sep = "")
ccss_org_ccss_exp.bim$KEY1 <- paste0(ccss_org_ccss_exp.bim$V1, ":", ccss_org_ccss_exp.bim$V4)
sum(duplicated(ccss_org_ccss_exp.bim$KEY1))
ccss_org_ccss_exp$REF <- ccss_org_ccss_exp.bim$V6[match(ccss_org_ccss_exp$SNP, ccss_org_ccss_exp.bim$V2)]
ccss_org_ccss_exp$ALT <- ccss_org_ccss_exp.bim$V5[match(ccss_org_ccss_exp$SNP, ccss_org_ccss_exp.bim$V2)]
table(ccss_org_ccss_exp$ALT == ccss_org_ccss_exp$A1)
ccss_org_ccss_exp[ccss_org_ccss_exp$ALT != ccss_org_ccss_exp$A1, c("ALT", "REF")] <- ccss_org_ccss_exp[ccss_org_ccss_exp$ALT != ccss_org_ccss_exp$A1, c("REF", "ALT")]
ccss_org_ccss_exp$OR_95CI <- paste0(round(ccss_org_ccss_exp$OR,2), " (", round(ccss_org_ccss_exp$L95, 2), "-", round(ccss_org_ccss_exp$U95,2), ")")
ccss_org_ccss_exp$KEY <- paste0("chr", ccss_org_ccss_exp$CHR, ":", ccss_org_ccss_exp$BP, ":", ccss_org_ccss_exp$REF, ":", ccss_org_ccss_exp$A1)  
ccss_org_ccss_exp <- ccss_org_ccss_exp[c("OR_95CI", "P", "KEY")]


library(purrr)

dfs <- list(combined, sjlife, ccss_org, ccss_exp, ccss_org_ccss_exp)
final_merged_df <- reduce(dfs, merge, by = "KEY")
colnames(final_merged_df) <- c("SNP", "KEY", "CHR", "BP", "A1", "REF", "OR_95CI_combined", "P_combined", 
                               "OR_95CI_sjlife", "P_sjlife", 
                               "OR_95CI_ccss_org", "P_ccss_org", 
                               "OR_95CI_ccss_exp", "P_ccss_exp", 
                               "OR_95CI_ccss_org_exp", "P_ccss_org_exp")

freq.df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples_updated_freq_out.frq", header = T)
dim(freq.df)
final_merged_df$combinedMAF <- freq.df$MAF[match(final_merged_df$SNP, freq.df$SNP)]
final_merged_df <- final_merged_df[final_merged_df$combinedMAF > 0.05,]

# read annoatation file
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

annotation <- annotation[grepl("missense", annotation$EFFECT, ignore.case = T),]
sum(final_merged_df$KEY %in% annotation$SNPID)
cc <- final_merged_df[final_merged_df$KEY %in% annotation$SNPID,]
sum(duplicated(cc$KEY))
cc$variation <- annotation$EFFECT[match(cc$KEY, annotation$SNPID)]
cc$Gene <- annotation$Gene[match(cc$KEY, annotation$SNPID)]


write.table(cc, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/common_missense_variants.txt", col.names = T, sep = "\t", quote = F, row.names = F)
