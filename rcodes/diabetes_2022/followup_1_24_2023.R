load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/Results_final_to_Cindy/SJLIFE_T2D_GWAS_data.RData")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/followup_jan_24_2023")

AFR <- dat_afr[1:2]
AFR$FID <- 0

EUR <- dat_eur[1:2]
EUR$FID <- 0

write.table(AFR, "AFR_samples.list", sep = " ", quote = FALSE, row.names = F)
write.table(EUR, "EUR_samples.list", sep = " ", quote = FALSE, row.names = F)



######################
## Get genotype AFR ##
######################
raw <- read.table("folowup_24_2023_AFR_recodeA.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[!grepl("FID|IID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
colnames(raw) = gsub("\\.", ":", HEADER)
raw$IID <- row.names(raw)

sum(raw$IID %in% dat_afr$IID)
# 574

raw <- cbind.data.frame(dat_afr[match(raw$IID, dat_afr$IID),3:23], raw)

# mtcars %>%
#   relocate(gear, carb, .before = cyl) %>%
#   head()

raw <- raw %>% 
  relocate(IID, .before = t2d)

AFR_six_vars <- raw

######################
## Get genotype EUR ##
######################
raw <- read.table("folowup_24_2023_EUR_recodeA.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[!grepl("FID|IID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
colnames(raw) = gsub("\\.", ":", HEADER)
raw$IID <- row.names(raw)

sum(raw$IID %in% dat_eur$IID)
# 3102

raw <- cbind.data.frame(dat_eur[match(raw$IID, dat_eur$IID),3:23], raw)

# mtcars %>%
#   relocate(gear, carb, .before = cyl) %>%
#   head()

raw <- raw %>% 
  relocate(IID, .before = t2d)

EUR_six_vars <- raw

save.image("followup_jan_24_2023.RData")

df <- read.table(text = ("Chr	BP	EA	NEA
8	50635641	T	C
8	50636999	C	T
8	50637285	A	T
8	50637444	A	G
8	50638825	G	A
8	50641552	A	G
8	50643118	T	G
8	50644980	A	G
8	50660409	G	GAAAA
8	50661918	G	C
8	50666979	A	G
8	50667110	A	G
8	50667850	G	T
8	50671413	C	T
8	50671938	A	G
8	50672410	T	A
8	50672830	C	T
8	50673224	G	A
8	50674706	G	A
8	50675283	G	A
8	50676103	C	T
8	50676751	G	T
8	50681567	G	A"), header = T)

df$SNP <- paste0("chr", df$Chr, ":", df$BP, ":",df$NEA, ":", df$EA)

metasoft.res <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/METASOFT/metasoft_res_edited", sep = "\t", header=F)
header <- c("RSID", "STUDY", "PVALUE_FE", "BETA_FE", "STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2")
# header <- c("RSID", "STUDY", "PVALUE_FE", "BETA_FE", "STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE", "PVALUES_OF_STUDIES", "MVALUES_OF_STUDIES")
length(header)
metasoft.res <- metasoft.res[,1:11]
colnames(metasoft.res) <- header 
head(metasoft.res)

sum(df$SNP %in% metasoft.res$RSID)

df <- cbind.data.frame(df, metasoft.res[match(df$SNP, metasoft.res$RSID),])
rm(metasoft.res)
write.table(df, "followup_jan24_2023_metasoft_res", sep = " ", quote = FALSE, row.names = F)
