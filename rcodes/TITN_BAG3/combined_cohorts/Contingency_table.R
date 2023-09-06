## Contingency_table
# 1. --------------------------- TTN
PHENO <- read.table("pheno/ccss_exp_eur_cardiotoxic_exposed.pheno", header = T)
# chr2_178562809_T_C Before removing EOAD samples
before <- read.table("ccss_exp_before_chr2_178562809_T_C.raw", header = T)
PHENO$geno.before <- before$chr2.178562809.T.C_C[match(PHENO$FID, before$FID)]

table(PHENO$CMP2plus, PHENO$geno.before)
#    0   1   2
# 1 876 481  67
# 2  34  22   4


# chr2_178562809_T_C After removing EOAD samples
after <- read.table("ccss_exp_after_chr2_178562809_T_C.raw", header = T)
PHENO$geno.after <- after$chr2.178562809.T.C_C[match(PHENO$FID, after$FID)]
PHENO <- PHENO[!is.na(PHENO$geno.after),]
table(PHENO$CMP2plus, PHENO$geno.after)
#    0   1   2
# 1 876 481  67
# 2  18  11   4

# 2. --------------------------- BAG3
# chr10_119670121_T_C
PHENO <- read.table("pheno/ccss_exp_eur_cardiotoxic_exposed.pheno", header = T)
# chr2_178562809_T_C Before removing EOAD samples
before <- read.table("ccss_exp_before_chr10_119670121_T_C.raw", header = T)
PHENO$geno.before <- before$chr10.119670121.T.C_C[match(PHENO$FID, before$FID)]

table(PHENO$CMP2plus, PHENO$geno.before)
# 0   1   2
# 1 879 481  64
# 2  45  13   2


# chr2_178562809_T_C After removing EOAD samples
after <- read.table("ccss_exp_after_chr10_119670121_T_C.raw", header = T)
PHENO$geno.after <- after$chr10.119670121.T.C_C[match(PHENO$FID, after$FID)]
PHENO <- PHENO[!is.na(PHENO$geno.after),]
table(PHENO$CMP2plus, PHENO$geno.after)
# 0   1   2
# 1 879 481  64
# 2  24   8   1

