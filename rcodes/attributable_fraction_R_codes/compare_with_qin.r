##################################
## Sanity check with Qin's data ##
##################################
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/cleaned_phenotype.RDATA")
pheno.crosscheck <- PHENO.ANY_SN
pheno.crosscheck <- pheno.crosscheck[mixedorder(pheno.crosscheck$sjlid),]
qi.df <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/toyadav.sas7bdat")
qi.df <- qi.df[mixedorder(qi.df$sjlid),]
head(qi.df)
qi.df.unique <- qi.df[!duplicated(qi.df$sjlid),]
# length(unique(qi.df.unique$sjlid))
# # 4402
# sum(pheno.crosscheck$sjlid %in% qi.df.unique$sjlid)
# # 4402
# qi.df <- qi.df[!duplicated(qi.df$sjlid),]
table(qi.df$MutationStatus)
table(qi.df$sn)
# 1263

qi.df.SN <- qi.df[grepl(1, qi.df$sn),]
dim(qi.df.SN)
# 1263
# Keeping unique rows
qi.df.SN <- distinct(qi.df.SN)
dim(qi.df.SN)



qi.df.not.SN <- qi.df[!grepl(1, qi.df$sn),]
dim(qi.df.not.SN)
# 3911
# 4402-3911 = 491

# Extracting SNs for the first time from Qi's data
qi.df.SN.filtered <- setDT(qi.df.SN)[,.SD[which.min(evaldt)],by=sjlid][order(evaldt, decreasing = FALSE)]
dim(qi.df.SN.filtered)

# Now merge back the data with SN and without SN 
qi.df.final <- rbind.data.frame(qi.df.not.SN, qi.df.SN.filtered)
dim(qi.df.final)
qi.df.final <- qi.df.final[mixedorder(qi.df.final$sjlid),]

paste0(table(qi.df.final$anyrt_5)[2], ", ", table(pheno.crosscheck$anyrt_5)[2]) # anyrt_5
# "2169, 2174"
paste0(table(qi.df.final$brainrt_yn)[2], ", ", table(pheno.crosscheck$brainrt_yn)[2]) # Brain
# "858, 1095"
paste0(table(qi.df.final$neckrt_yn)[2], ", ", table(pheno.crosscheck$neckrt_yn)[2]) # Neck
# "601, 797"
paste0(table(qi.df.final$chest_yn)[2], ", ", table(pheno.crosscheck$chestrt_yn)[2]) # Chest
# 628, 861
paste0(table(qi.df.final$abdomenrt_yn)[2], ", ", table(pheno.crosscheck$abdomenrt_yn)[2]) # Abdomen
# "567, 790"
paste0(table(qi.df.final$pelvis_yn)[2], ", ", table(pheno.crosscheck$pelvisrt_yn)[2]) # Pelvis
# "506, 688"
paste0(table(qi.df.final$aaclassic_5 ==1)[2], ", ", table(pheno.crosscheck$aa_class_dose_5_yn)[2]) # Alkylating_5
# "2480, 2473"


paste0(table(qi.df.final$anthracyclines_5 > 0)[2], ", ", table(pheno.crosscheck$anthra_jco_dose_5_yn)[2]) # Anthracyclines_5
# "2455, 2451"
paste0(table(qi.df.final$epipodophyllotoxins_5 > 0)[2], ", ", table(pheno.crosscheck$epitxn_dose_5_yn)[2]) # Epipodophyllotoxins_5
# "1496, 1488"

# Age at diagnosis
paste0(paste0(paste0( median(qi.df.final$agedx),  " (", paste0( quantile(qi.df.final$agedx, 1 / 4), "−",
                quantile(qi.df.final$agedx, 3 / 4)), ")" )), ", ", paste0(paste0(
                round(median(pheno.crosscheck$agedx), 2), " (", paste0(round(quantile( pheno.crosscheck$agedx, 1 / 4 ), 2), "−",
                round(quantile( pheno.crosscheck$agedx, 3 / 4), 2), ")")))) 
# "6.275 (2.79−12.525), 6.28 (2.79−12.52)"


table(qi.df.final$diaggrp == pheno.crosscheck$diaggrp)
# FALSE  TRUE 
# 1087  3315

table(qi.df.final$epitxn_dose)
sum(is.na(qi.df.final$epipodophyllotoxins_5))
# 10 
qi.df.final[is.na(qi.df.final$epipodophyllotoxins_5),c("sjlid", "epipodophyllotoxins_5")]

pheno.crosscheck.SN.491 <- pheno.crosscheck[pheno.crosscheck$sjlid %in%  qi.df.SN.filtered$sjlid,]
pheno.crosscheck.SN.491 <- pheno.crosscheck.SN.491[mixedorder(pheno.crosscheck.SN.491$sjlid),]
table(pheno.crosscheck.SN.491$diaggrp == qi.df.SN.filtered$diaggrp)
# FALSE  TRUE 
# 376   115 

pheno.crosscheck.not.SN <- pheno.crosscheck[pheno.crosscheck$sjlid %in%  qi.df.not.SN$sjlid,]
pheno.crosscheck.not.SN <- pheno.crosscheck.not.SN[mixedorder(pheno.crosscheck.not.SN$sjlid),]
table(pheno.crosscheck.not.SN$diaggrp == qi.df.not.SN$diaggrp)