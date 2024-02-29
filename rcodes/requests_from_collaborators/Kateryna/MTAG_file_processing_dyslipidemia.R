library(data.table)
#################
## Process EUR ##
#################

## 1. chltot
eur_chltot <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_eur_chltot_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
eur_chltot.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_eur_chltot_frequency_edited.frq", header =T)
sum(eur_chltot.frq$MAF < 0.01, na.rm = T)
# 18278
## exclude rare
eur_chltot.frq <- eur_chltot.frq[which(eur_chltot.frq$MAF >= 0.01),]

eur_chltot$MAF <- eur_chltot.frq$MAF[match(eur_chltot$SNP, eur_chltot.frq$SNP)]
eur_chltot <- eur_chltot[!is.na(eur_chltot$MAF),]
# > head(eur_chltot)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

eur_chltot <- eur_chltot[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(eur_chltot) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

## 2. hdl
eur_hdl <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_eur_hdl_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
eur_hdl.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_eur_hdl_frequency_edited.frq", header =T)
sum(eur_hdl.frq$MAF < 0.01, na.rm = T)
# 20457
## exclude rare
eur_hdl.frq <- eur_hdl.frq[which(eur_hdl.frq$MAF >= 0.01),]

eur_hdl$MAF <- eur_hdl.frq$MAF[match(eur_hdl$SNP, eur_hdl.frq$SNP)]
eur_hdl <- eur_hdl[!is.na(eur_hdl$MAF),]
# > head(eur_hdl)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

eur_hdl <- eur_hdl[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(eur_hdl) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

## 3. ldl
eur_ldl <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_eur_ldl_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
eur_ldl.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_eur_ldl_frequency_edited.frq", header =T)
sum(eur_ldl.frq$MAF < 0.01, na.rm = T)
# 19465
## exclude rare
eur_ldl.frq <- eur_ldl.frq[which(eur_ldl.frq$MAF >= 0.01),]

eur_ldl$MAF <- eur_ldl.frq$MAF[match(eur_ldl$SNP, eur_ldl.frq$SNP)]
eur_ldl <- eur_ldl[!is.na(eur_ldl$MAF),]
# > head(eur_ldl)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

eur_ldl <- eur_ldl[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(eur_ldl) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

## 4. nonhdl
eur_nonhdl <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_eur_nonhdl_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
eur_nonhdl.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_eur_nonhdl_frequency_edited.frq", header =T)
sum(eur_nonhdl.frq$MAF < 0.01, na.rm = T)
# 19533
## exclude rare
eur_nonhdl.frq <- eur_nonhdl.frq[which(eur_nonhdl.frq$MAF >= 0.01),]

eur_nonhdl$MAF <- eur_nonhdl.frq$MAF[match(eur_nonhdl$SNP, eur_nonhdl.frq$SNP)]
eur_nonhdl <- eur_nonhdl[!is.na(eur_nonhdl$MAF),]
# > head(eur_nonhdl)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

eur_nonhdl <- eur_nonhdl[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(eur_nonhdl) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

## 5. trigly
eur_trigly <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_eur_trigly_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
eur_trigly.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_eur_trigly_frequency_edited.frq", header =T)
sum(eur_trigly.frq$MAF < 0.01, na.rm = T)
# 20140
## exclude rare
eur_trigly.frq <- eur_trigly.frq[which(eur_trigly.frq$MAF >= 0.01),]

eur_trigly$MAF <- eur_trigly.frq$MAF[match(eur_trigly$SNP, eur_trigly.frq$SNP)]
eur_trigly <- eur_trigly[!is.na(eur_trigly$MAF),]
# > head(eur_trigly)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

eur_trigly <- eur_trigly[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(eur_trigly) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")


# We have dataframes named eur_chltot, eur_ldl, eur_hdl, eur_trigly, eur_nonhdl
# Create a vector with the overlapping variants 
common_variants <- Reduce(intersect, list(eur_chltot$snpid, eur_ldl$snpid, eur_hdl$snpid, eur_trigly$snpid, eur_nonhdl$snpid))

get.rsid <- eur_chltot[eur_chltot$snpid %in% common_variants,]
get.rsid$fullSNP <- paste0("chr",get.rsid$snpid, ":", get.rsid$a2, ":",get.rsid$a1)

# Filter each dataframe to keep only the common variants
eur_chltot.shared <- eur_chltot[eur_chltot$snpid %in% common_variants, ]
eur_ldl.shared <- eur_ldl[eur_ldl$snpid %in% common_variants, ]
eur_hdl.shared <- eur_hdl[eur_hdl$snpid %in% common_variants, ]
eur_trigly.shared <- eur_trigly[eur_trigly$snpid %in% common_variants, ]
eur_nonhdl.shared <- eur_nonhdl[eur_nonhdl$snpid %in% common_variants, ]

write.table(eur_chltot.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/eur_chltot.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(eur_ldl.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/eur_ldl.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(eur_hdl.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/eur_hdl.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(eur_trigly.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/eur_trigly.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(eur_nonhdl.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/eur_nonhdl.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#################
## Process AFR ##
#################

## 1. chltot
afr_chltot <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_afr_chltot_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
afr_chltot.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_afr_chltot_frequency_edited.frq", header =T)
sum(afr_chltot.frq$MAF < 0.01, na.rm = T)
# 13985
## exclude rare
afr_chltot.frq <- afr_chltot.frq[which(afr_chltot.frq$MAF >= 0.01),]

afr_chltot$MAF <- afr_chltot.frq$MAF[match(afr_chltot$SNP, afr_chltot.frq$SNP)]
afr_chltot <- afr_chltot[!is.na(afr_chltot$MAF),]
# > head(afr_chltot)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

afr_chltot <- afr_chltot[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(afr_chltot) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

## 2. hdl
afr_hdl <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_afr_hdl_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
afr_hdl.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_afr_hdl_frequency_edited.frq", header =T)
sum(afr_hdl.frq$MAF < 0.01, na.rm = T)
# 15614
## exclude rare
afr_hdl.frq <- afr_hdl.frq[which(afr_hdl.frq$MAF >= 0.01),]

afr_hdl$MAF <- afr_hdl.frq$MAF[match(afr_hdl$SNP, afr_hdl.frq$SNP)]
afr_hdl <- afr_hdl[!is.na(afr_hdl$MAF),]
# > head(afr_hdl)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

afr_hdl <- afr_hdl[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(afr_hdl) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

## 3. ldl
afr_ldl <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_afr_ldl_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
afr_ldl.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_afr_ldl_frequency_edited.frq", header =T)
sum(afr_ldl.frq$MAF < 0.01, na.rm = T)
# 14813
## exclude rare
afr_ldl.frq <- afr_ldl.frq[which(afr_ldl.frq$MAF >= 0.01),]

afr_ldl$MAF <- afr_ldl.frq$MAF[match(afr_ldl$SNP, afr_ldl.frq$SNP)]
afr_ldl <- afr_ldl[!is.na(afr_ldl$MAF),]
# > head(afr_ldl)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

afr_ldl <- afr_ldl[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(afr_ldl) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

## 4. nonhdl
afr_nonhdl <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_afr_nonhdl_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
afr_nonhdl.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_afr_nonhdl_frequency_edited.frq", header =T)
sum(afr_nonhdl.frq$MAF < 0.01, na.rm = T)
# 14840
## exclude rare
afr_nonhdl.frq <- afr_nonhdl.frq[which(afr_nonhdl.frq$MAF >= 0.01),]

afr_nonhdl$MAF <- afr_nonhdl.frq$MAF[match(afr_nonhdl$SNP, afr_nonhdl.frq$SNP)]
afr_nonhdl <- afr_nonhdl[!is.na(afr_nonhdl$MAF),]
# > head(afr_nonhdl)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

afr_nonhdl <- afr_nonhdl[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(afr_nonhdl) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

## 5. trigly
afr_trigly <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/DYSLPDM_afr_trigly_chrALL.assoc.linear.clean.Psorted.formetal.common", header =T)
#read frequency file
afr_trigly.frq <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/freq/DYSLPDM_afr_trigly_frequency_edited.frq", header =T)
sum(afr_trigly.frq$MAF < 0.01, na.rm = T)
# 15523
## exclude rare
afr_trigly.frq <- afr_trigly.frq[which(afr_trigly.frq$MAF >= 0.01),]

afr_trigly$MAF <- afr_trigly.frq$MAF[match(afr_trigly$SNP, afr_trigly.frq$SNP)]
afr_trigly <- afr_trigly[!is.na(afr_trigly$MAF),]
# > head(afr_trigly)
#     CHR         SNP       BP A1 TEST NMISS       OR      SE     L95     U95   STAT         P REF ALT A2      BETA     MAF
# 1:  19 19:44908822 44908822  T  ADD   497 0.568303 0.09024 -0.7420 -0.3882 -6.262 8.420e-10   C   T  C -0.565101 0.07492
# 2:  19 19:44897490 44897490  A  ADD   498 0.592977 0.09166 -0.7022 -0.3429 -5.701 2.077e-08   T   A  T -0.522600 0.03120

afr_trigly <- afr_trigly[,c("SNP", "CHR", "BP", "A1", "A2", "MAF", "STAT", "P", "NMISS")]
colnames(afr_trigly) <- c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")

# We have dataframes named afr_chltot, afr_ldl, afr_hdl, afr_trigly, afr_nonhdl
# Create a vector with the overlapping variants 
common_variants <- Reduce(intersect, list(afr_chltot$snpid, afr_ldl$snpid, afr_hdl$snpid, afr_trigly$snpid, afr_nonhdl$snpid))

get.rsid <- afr_chltot[afr_chltot$snpid %in% common_variants,]
get.rsid$fullSNP <- paste0("chr",get.rsid$snpid, ":", get.rsid$a2, ":",get.rsid$a1)

# Filter each dataframe to keep only the common variants
afr_chltot.shared <- afr_chltot[afr_chltot$snpid %in% common_variants, ]
afr_ldl.shared <- afr_ldl[afr_ldl$snpid %in% common_variants, ]
afr_hdl.shared <- afr_hdl[afr_hdl$snpid %in% common_variants, ]
afr_trigly.shared <- afr_trigly[afr_trigly$snpid %in% common_variants, ]
afr_nonhdl.shared <- afr_nonhdl[afr_nonhdl$snpid %in% common_variants, ]

write.table(afr_chltot.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/afr_chltot.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(afr_ldl.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/afr_ldl.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(afr_hdl.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/afr_hdl.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(afr_trigly.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/afr_trigly.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(afr_nonhdl.shared, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Dyslipidemia/GWAS/MTAG/afr_nonhdl.shared.txt", col.names = T, row.names = F, sep = "\t", quote = F)

