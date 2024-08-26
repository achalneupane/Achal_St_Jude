############################################################
## Prepare file to run PRSice with shell script PRSice.sh ##
############################################################
# steps:
# 1. Download PRS files from https://www.pgscatalog.org/ for GRCh38 build. For example, PGS000688 can be downloaded from https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000688/ScoringFiles/Harmonized/
# 2. Prepare PRS files so that it has CHR, GRCh38 BP (Base position), SNP ID, other_allele (Ref), effect_allele (Alt), BETA and P. If P is missing, you can set P=1.
# 3. Once the *hmPOS_GRCh38_processed.txt file is prepared, proceed to run PRSice.sh
  
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_exp")
files <- list.files(pattern = "PGS.*hmPOS_GRCh38.txt$") ## Look for all Prs files in the folder that need to be processed.

# Here, I am reading each PRS file and writing the processed output file in loop in ccss_exp folder
for (i in 1:length(files)){
print(paste0("Doing: ", files[i]))
df <-   read.table(files[i], sep = "\t", header = T)
# P <- runif(nrow(df), min = 0, max = 1)
P <- 1
cc <- cbind.data.frame(CHR=df$chr_name, BP=df$hm_pos, 
                       SNP=gsub(" ","", paste0("chr",df$chr_name,":", df$hm_pos, ":", df$other_allele, ":", df$effect_allele)), 
                       other_allele=df$other_allele, effect_allele=df$effect_allele, BETA=df$effect_weight, P=P)
cc <- cc[(!duplicated(cc$SNP)),]
cc <- cc[!(is.na(cc$BP)),]
writefile <- gsub(".txt","_processed.txt", files[i])
write.table(cc, writefile, sep = "\t", col.names = T, row.names = F, quote = F)
}


setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_org")
files <- list.files(pattern = "PGS.*hmPOS_GRCh38.txt$")

# Here, I am reading each PRS file and writing the processed output file in loop in ccss_org folder
for (i in 1:length(files)){
  print(paste0("Doing: ", files[i]))
  df <-   read.table(files[i], sep = "\t", header = T)
  # P <- runif(nrow(df), min = 0, max = 1)
  P <- 1
  cc <- cbind.data.frame(CHR=df$chr_name, BP=df$chr_position, 
                         SNP=gsub(" ","", paste0("chr",df$chr_name,":", df$chr_position, ":", df$other_allele, ":", df$effect_allele)), 
                         other_allele=df$other_allele, effect_allele=df$effect_allele, BETA=df$effect_weight, P=P)
  cc <- cc[(!duplicated(cc$SNP)),]
  cc <- cc[!(is.na(cc$BP)),]
  writefile <- gsub(".txt","_processed.txt", files[i])
  write.table(cc, writefile, sep = "\t", col.names = T, row.names = F, quote = F)
}

## Now run PRSice.sh


# # This is an additional comparison I was conducting to verify if the results align with the Plink PRS calculations, so there's no need to run it.
# ####################################################################
# ## Check whether PRSICE and PLINK methods give correlated results ##
# ####################################################################
# plink.method <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_exp/prs_out/PGS000688_prs_tab_separated.profile", header = T, sep = "\t")
# prsice.method <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/ccss_exp/PGS000688_PRS_output.all_score", header = T, sep = " ")
# prsice.method$plinkscores <- plink.method$SCORE[match(prsice.method$IID, plink.method$IID)]
# # Compute the correlation
# correlation <- cor(prsice.method$Pt_1, prsice.method$plinkscores)
# correlation
# # 0.9955616
# print(paste("Correlation: ", correlation))
# 
# # Create a scatter plot with a regression line
# library(ggplot2)
# ggplot(prsice.method, aes(x = Pt_1, y = plinkscores)) +
#   geom_point() +
#   geom_smooth(method = "lm", col = "blue") +
#   labs(title = paste("Correlation between Pt_1 and Plink Scores: ", round(correlation, 2)),
#        x = "Pt_1",
#        y = "Plink Scores") +
#   theme_minimal()



