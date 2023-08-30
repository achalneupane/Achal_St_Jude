## Read association results
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3")
df.asso <- read.table("sjlife_ccss_org_ccss_exp_samples_results_v2.assoc.logistic", header = T, stringsAsFactors = F)
dim(df.asso)

df.asso$bim_SNP <- df.asso$SNP 
df.asso$SNP <- sapply(strsplit(df.asso$SNP, split=";"), "[", 1)

df.asso$tmp.A2 <- gsub("DEL", ":DEL", sapply(strsplit(gsub(":DEL", "DEL", df.asso$SNP), split=":"), "[", 3))
df.asso$tmp.A1 <- gsub("DEL", ":DEL", sapply(strsplit(gsub(":DEL", "DEL", df.asso$SNP), split=":"), "[", 4))


df.asso$A2 <- NA
df.asso$A2[df.asso$A1 == df.asso$tmp.A1] <- df.asso$tmp.A2[df.asso$A1 == df.asso$tmp.A1]
df.asso$A2[df.asso$A1 != df.asso$tmp.A1] <- df.asso$tmp.A1[df.asso$A1 != df.asso$tmp.A1]

df.asso$KEY <- paste0("chr", df.asso$CHR, ":", df.asso$BP) 
df.asso$SNP <- paste0("chr", df.asso$CHR, ":", df.asso$BP, ":", df.asso$A2, ":", df.asso$A1) 

df.asso$REF <- df.asso$A2
df.asso$ALT <- df.asso$A1
## Flip alleles
df.asso$REF_flipped <- chartr("acgtACGT", "tgcaTGCA", df.asso$REF)
df.asso$ALT_flipped <- chartr("acgtACGT", "tgcaTGCA", df.asso$ALT)

BAG3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/chr10_119_list.txt", header = T)
TTN <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/chr2_178_list.txt", header = T)
TTN_BAG3 <- rbind.data.frame(BAG3,TTN)
TTN_BAG3$KEY <- paste(TTN_BAG3$CHROM, TTN_BAG3$POS, sep = ":")
TTN_BAG3$SNP <- paste(TTN_BAG3$CHROM, TTN_BAG3$POS, TTN_BAG3$REF, TTN_BAG3$ALT, sep = ":")

## KEEP missense only
TTN_BAG3 <- TTN_BAG3[grepl("missense", TTN_BAG3$ANN....EFFECT, ignore.case = T),]
## Remove duplicates
TTN_BAG3 <- TTN_BAG3[!duplicated(TTN_BAG3$SNP),]
dim(TTN_BAG3)
# 3007    25

# # TTN_BAG3$SNP <- sub(";.*", "", TTN_BAG3$ID)
# sum(df.asso$KEY %in% TTN_BAG3$KEY)
## 111
# sum(TTN_BAG3$KEY %in% df.asso$KEY)
## 112
sum(df.asso$SNP %in% TTN_BAG3$SNP)
# 104




for (i in 1:nrow(df.asso)){
  print(paste0("Doing iteration: ", i))
  if (sum(df.asso$KEY[i] %in% TTN_BAG3$KEY) > 0){ # Only if position matches; do
    match.index <- grep(df.asso$KEY[i], TTN_BAG3$SNP)
    for(j in 1:length(match.index)){ # direct match or match by swapping alleles
      if(df.asso$REF[i] == TTN_BAG3$REF[match.index[j]] & df.asso$ALT[i] == TTN_BAG3$ALT[match.index[j]]){
         df.asso$MATCH_gnomAD[i] <- "DIRECT_MATCH"
         df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else if
        (df.asso$REF[i] == TTN_BAG3$ALT[match.index[j]] & df.asso$ALT[i] == TTN_BAG3$REF[match.index[j]]){
        df.asso$MATCH_gnomAD[i] <- "INDIRECT_MATCH0"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else if # match by flipping alleles
      (df.asso$REF_flipped[i] == TTN_BAG3$REF[match.index[j]] & df.asso$ALT_flipped[i] == TTN_BAG3$ALT[match.index[j]]){
        df.asso$MATCH_gnomAD[i] <- "INDIRECT_MATCH1"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else if # match by swapping flipped alleles
      (df.asso$REF_flipped[i] == TTN_BAG3$ALT[match.index[j]] & df.asso$ALT_flipped[i] == TTN_BAG3$REF[match.index[j]]){
        df.asso$MATCH_gnomAD[i] <- "INDIRECT_MATCH2"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else if # match by flipping one allele or by swapping one of the flipped alleles
      ((df.asso$REF_flipped[i] == TTN_BAG3$REF[match.index[j]] & df.asso$ALT[i] == TTN_BAG3$ALT[match.index[j]])|
       (df.asso$ALT_flipped[i] == TTN_BAG3$REF[match.index[j]] & df.asso$REF[i] == TTN_BAG3$ALT[match.index[j]])){
        df.asso$MATCH_gnomAD[i] <- "INDIRECT_MATCH3"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else{
        # df.asso$MATCH_gnomAD[i] <- NA
        # df.asso$gnomAD_equivalent[i] <- NA
        # next
        print(i,j)
      }
      
    }
  } else {
    df.asso$MATCH_gnomAD[i] <- NA
    df.asso$gnomAD_equivalent[i] <- NA
  }
  
}


missense <- df.asso[!is.na(df.asso$MATCH_gnomAD),]
dim(missense)
## 110
missense$Substitution <- TTN_BAG3$ANN....EFFECT[match(missense$gnomAD_equivalent, TTN_BAG3$SNP)]
missense$Impact <- TTN_BAG3$ANN....IMPACT[match(missense$gnomAD_equivalent, TTN_BAG3$SNP)]
missense$Gene <- TTN_BAG3$ANN....GENE[match(missense$gnomAD_equivalent, TTN_BAG3$SNP)]
###############################################
## Now repeat this for non-missense variants ##
###############################################
df.asso <- df.asso[!df.asso$SNP %in% missense$SNP,]

BAG3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/chr10_119_list.txt", header = T)
TTN <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/chr2_178_list.txt", header = T)
TTN_BAG3 <- rbind.data.frame(BAG3,TTN)
TTN_BAG3$KEY <- paste(TTN_BAG3$CHROM, TTN_BAG3$POS, sep = ":")
TTN_BAG3$SNP <- paste(TTN_BAG3$CHROM, TTN_BAG3$POS, TTN_BAG3$REF, TTN_BAG3$ALT, sep = ":")

## KEEP nonmissense only
TTN_BAG3 <- TTN_BAG3[!grepl("missense", TTN_BAG3$ANN....EFFECT, ignore.case = T),]
dim(TTN_BAG3)
# 688003     25
## Remove duplicates
TTN_BAG3 <- TTN_BAG3[!duplicated(TTN_BAG3$SNP),]
dim(TTN_BAG3)
# 101732     25

# # TTN_BAG3$SNP <- sub(";.*", "", TTN_BAG3$ID)
# sum(df.asso$KEY %in% TTN_BAG3$KEY)
## 3386
# sum(TTN_BAG3$KEY %in% df.asso$KEY)
## 4140
sum(df.asso$SNP %in% TTN_BAG3$SNP)
# 2759




for (i in 1:nrow(df.asso)){
  print(paste0("Doing iteration: ", i))
  if (sum(df.asso$KEY[i] %in% TTN_BAG3$KEY) > 0){ # Only if position matches; do
    match.index <- grep(df.asso$KEY[i], TTN_BAG3$SNP)
    for(j in 1:length(match.index)){ # direct match or match by swapping alleles
      if(df.asso$REF[i] == TTN_BAG3$REF[match.index[j]] & df.asso$ALT[i] == TTN_BAG3$ALT[match.index[j]]){
        df.asso$MATCH_gnomAD[i] <- "DIRECT_MATCH"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else if
      (df.asso$REF[i] == TTN_BAG3$ALT[match.index[j]] & df.asso$ALT[i] == TTN_BAG3$REF[match.index[j]]){
        df.asso$MATCH_gnomAD[i] <- "INDIRECT_MATCH0"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else if # match by flipping alleles
      (df.asso$REF_flipped[i] == TTN_BAG3$REF[match.index[j]] & df.asso$ALT_flipped[i] == TTN_BAG3$ALT[match.index[j]]){
        df.asso$MATCH_gnomAD[i] <- "INDIRECT_MATCH1"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else if # match by swapping flipped alleles
      (df.asso$REF_flipped[i] == TTN_BAG3$ALT[match.index[j]] & df.asso$ALT_flipped[i] == TTN_BAG3$REF[match.index[j]]){
        df.asso$MATCH_gnomAD[i] <- "INDIRECT_MATCH2"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else if # match by flipping one allele or by swapping one of the flipped alleles
      ((df.asso$REF_flipped[i] == TTN_BAG3$REF[match.index[j]] & df.asso$ALT[i] == TTN_BAG3$ALT[match.index[j]])|
       (df.asso$ALT_flipped[i] == TTN_BAG3$REF[match.index[j]] & df.asso$REF[i] == TTN_BAG3$ALT[match.index[j]])){
        df.asso$MATCH_gnomAD[i] <- "INDIRECT_MATCH3"
        df.asso$gnomAD_equivalent[i] <- TTN_BAG3$SNP[match.index[j]]
      } else{
        # df.asso$MATCH_gnomAD[i] <- NA
        # df.asso$gnomAD_equivalent[i] <- NA
        # next
        print(i,j)
      }
      
    }
  } else {
    df.asso$MATCH_gnomAD[i] <- NA
    df.asso$gnomAD_equivalent[i] <- NA
  }
  
}


nonmissense <- df.asso[!is.na(df.asso$MATCH_gnomAD),]
nonmissense$Substitution <- TTN_BAG3$ANN....EFFECT[match(nonmissense$gnomAD_equivalent, TTN_BAG3$SNP)]
nonmissense$Impact <- TTN_BAG3$ANN....IMPACT[match(nonmissense$gnomAD_equivalent, TTN_BAG3$SNP)]
nonmissense$Gene <- TTN_BAG3$ANN....GENE[match(nonmissense$gnomAD_equivalent, TTN_BAG3$SNP)]
dim(nonmissense)
# 3051
ASSOC.results <- rbind.data.frame(missense, nonmissense)


## Add MAF
freq.df <- read.table("sjlife_ccss_org_ccss_exp_samples_freq_out.frq", header = T)
dim(freq.df)

sum(ASSOC.results$bim_SNP %in% freq.df$SNP)
# 3161
ASSOC.results$MAF_SJLIFE <- freq.df$MAF[match(ASSOC.results$bim_SNP, freq.df$SNP)]

ASSOC.results$new_OR <- paste0(ASSOC.results$OR, " (", ASSOC.results$L95, "-", ASSOC.results$U95, ")")


# write.table(ASSOC.results, "sjlife_ccss_org_ccss_exp__ttn_bag3.assoc.final_v2_07_06_2023.txt", col.names = T, row.names = F, quote = F, sep = "\t")

ASSOC.results <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp__ttn_bag3.assoc.final_v2_07_06_2023.txt", header = T, sep = "\t")

## ADD individual cohorts
SJLIFE <- read.table("sjlife_results_AN.assoc.logistic", header = T)
SJLIFE$new_OR <- paste0(SJLIFE$OR, " (", SJLIFE$L95, "-", SJLIFE$U95, ")")
SJLIFE$KEY <- paste0("chr", SJLIFE$CHR, ":", SJLIFE$BP)
SJLIFE$dups <- duplicated(SJLIFE$KEY)

CCSS_ORG <- read.table("ccss_org_results_AN.assoc.logistic", header = T)
CCSS_ORG$KEY <- paste0("chr", CCSS_ORG$CHR, ":", CCSS_ORG$BP)
CCSS_ORG$new_OR <- paste0(CCSS_ORG$OR, " (", CCSS_ORG$L95, "-", CCSS_ORG$U95, ")")

CCSS_EXP <- read.table("ccss_exp_results_AN.assoc.logistic", header = T)
CCSS_EXP$KEY <- paste0("chr", CCSS_EXP$CHR, ":", CCSS_EXP$BP)
CCSS_EXP$new_OR <- paste0(CCSS_EXP$OR, " (", CCSS_EXP$L95, "-", CCSS_EXP$U95, ")")

sum(!SJLIFE$KEY %in% ASSOC.results$KEY)

# These variants : rs3829746 and rs2234962 are not significant in CCSS, but only in SJLIFE, so decided to only include combined data
####























##################
## Forest plots ##
##################
df <- read.table(text = "SNP	rsID
chr2:178461008:C:G	rs17304212
chr2:178633315:T:C	rs6723526
chr2:178580212:C:T	rs2303838
chr2:178566270:G:A	rs3731746
chr2:178586693:G:A	rs2042996
chr2:178556967:A:G	rs9808377
chr2:178599800:T:C	rs1001238
chr2:178571293:G:A	rs744426
chr2:178532834:C:T	rs3829747
chr2:178562809:T:C	rs3829746
chr2:178541464:C:T	rs3731749
chr2:178592420:G:A	rs16866406
chr10:119670121:T:C	rs2234962
chr2:178593864:C:T	rs2288569
chr2:178693639:T:C	rs2042995", header = T)

sum(ASSOC.results$SNP %in% df$SNP)
ASSOC.results <- ASSOC.results[ASSOC.results$SNP %in% df$SNP,]
ASSOC.results$rsID <- df$rsID[match(ASSOC.results$SNP, df$SNP)]

ASSOC.results$lower <- ASSOC.results$L95
ASSOC.results$upper <- ASSOC.results$U95

summary_table <- ASSOC.results[c("rsID", "OR", "lower", "upper", "P")]

library(dplyr)
summary_table <- summary_table %>% arrange(factor(rsID, levels = df$rsID))

# library(ggplot2)
# P <- ggplot(data = summary_table, aes(x = OR, y = reorder(rsID, -P))) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   scale_x_continuous(limits = c(0.6, 2), breaks = seq(0.6, 2, 0.2)) +
#   xlab("OR (95% CI)") +
#   ylab("rsID") +
#   scale_y_discrete(labels = c("rs17304212", "rs6723526", "rs2303838", expression(bold(rs3829746)), "rs2042996", "rs9808377",
#                               "rs1001238", "rs744426", "rs3829747", "rs3829746", "rs3731749", "rs16866406",
#                               expression(bold(rs2234962)), "rs2288569", "rs2042995")) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank())


library(gridExtra)

custom_face <- ifelse(summary_table$rsID %in% c("rs3829746", "rs2234962"), "bold", "plain")

P <- ggplot(data = summary_table, aes(x = OR, y = reorder(rsID, -P))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0.6, 2), breaks = seq(0.6, 2, 0.2)) +
  xlab("OR (95% CI)") +
  ylab("") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(face = rev(custom_face), color = "black")
  )

P
ggsave("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Figures/common_variant_analysis.tiff", P, dpi = 600, width = 5, height = 3, units = "in")
