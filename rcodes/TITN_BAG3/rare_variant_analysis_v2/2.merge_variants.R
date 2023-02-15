df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Summary_results.txt", sep = "\t", header = T)
dim(df)


# read metal
# df.metal <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Meta_analysis_sjlife_ccss_org_ccss_exp_fixed_1.tbl", header = T)
df.metal <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Meta_analysis_sjlife_ccss_org_ccss_exp_after_5yrs_elapsed_age_fixed_1.tbl", header = T)
dim(df.metal)


head(df)
head(df.metal)
df.metal$Allele1 <- toupper(df.metal$Allele1)
df.metal$Allele2 <- toupper(df.metal$Allele2)

sum(df$SNP %in% df.metal$MarkerName) # All match
new.df <- cbind.data.frame(df,df.metal[match(df$SNP, df.metal$MarkerName), c("Allele1", "Allele2",  "Effect", "StdErr", "P.value")])

# View(new.df)
new.df$new_metaleffect <- NA
new.df$new_metaleffect[new.df$A1 == new.df$Allele1] <- new.df$Effect[new.df$A1 == new.df$Allele1]
new.df$new_metaleffect[new.df$A1 != new.df$Allele1] <- -new.df$Effect[new.df$A1 != new.df$Allele1]
# format(round(exp(new.df$Effect[new.df$A1 == new.df$Allele1]),2), nsmall = 2)

new.df$new_metal.OR <- round(exp(new.df$new_metaleffect), 2)

## lower CI of beta = beta - (se) x 1.96
# new.df$LCI.metal.OR <- format(round(exp(new.df$new_metaleffect - new.df$StdErr* 1.96), 2), nsmall = 2)
# new.df$UCI.metal.OR <- format(round(exp(new.df$new_metaleffect + new.df$StdErr* 1.96), 2), nsmall = 2)

new.df$LCI.metal.OR <- round(exp(new.df$new_metaleffect - new.df$StdErr* 1.96), 2)
new.df$UCI.metal.OR <- round(exp(new.df$new_metaleffect + new.df$StdErr* 1.96), 2)

new.df$new_metal.OR <- paste0(new.df$new_metal.OR, " (", new.df$LCI.metal.OR, "-", new.df$UCI.metal.OR, ")" )
FINAL <- new.df
## Now add New results from CCSS Exp
df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ccss_exp_results_with_elapsed_age_5yrs.assoc.logistic", header = T)
dim(df)
df$SNP1 <- paste(df$CHR, df$BP,df$A1, sep = ":")
new.df$SNP1 <- paste(new.df$SNP, new.df$A1, sep = ":")
new.df$SNP2 <- paste(new.df$SNP, new.df$A2, sep = ":")
sum(new.df$SNP1 %in% df$SNP1)
sum(new.df$SNP2 %in% df$SNP1)
df.a <- df[df$SNP1 %in% new.df$SNP1,]
df.b <- df[df$SNP1 %in% new.df$SNP2,]
df <- rbind.data.frame(df.a, df.b)

df$new_OR_ccss_exp <- paste0(df$OR, " (", df$L95, "-", df$U95, ")")
df$new_ccss_exp_P <- df$P
df$VARID <- paste(df$CHR, df$BP, sep = ":")

sum(FINAL$SNP %in% df$VARID)

FINAL$new_OR_ccss_exp <- df$new_OR_ccss_exp[match(FINAL$SNP, df$VARID)]
FINAL$new_ccss_exp_P <- df$new_ccss_exp_P[match(FINAL$SNP, df$VARID)]
write.table(FINAL, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/merged_new_results_from_Meta_analysis_sjlife_ccss_org_ccss_exp_after_5yrs_elapsed_age_fixed_1.txt", col.names = T, sep = "\t", quote = F, row.names = F)
