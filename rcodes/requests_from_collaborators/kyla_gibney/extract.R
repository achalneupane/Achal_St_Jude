df <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/request/Kyla/Within_50m.txt", header = T, sep = "\t")
dim(df)
head(df)


# read raw 
mapfile <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt", header = F, sep = " ", stringsAsFactors = F)
raw <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS/plink/kyla_gibney/kyla_gibney_vars.raw", sep = " ", header = T, check.names = F)
dim(raw)
head(raw)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

raw$rs6971 <- raw$`chr22:43162920:A:G`
raw$rs6971[raw$rs6971 == 0 ] <- "AA"
raw$rs6971[raw$rs6971 == 1 ] <- "AG"
raw$rs6971[raw$rs6971 == 2 ] <- "GG"
table(raw$rs6971)

df$WGS_ID <- mapfile$V2[match(df$SJLID, mapfile$V4)]
df$WGS_ID[df$SJLID == "SJL5144399"] <- "SJNORM045714_G1-TB-16-03948"

df$rs6971 <- raw$rs6971[match(df$WGS_ID, raw$IID)]



raw <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/kyla_gibney/kyla_gibney_vars.raw", sep = " ", header = T, check.names = F)
dim(raw)
head(raw)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

raw$rs6971 <- raw$`chr22:43162920:A:G`
raw$rs6971[raw$rs6971 == 0 ] <- "AA"
raw$rs6971[raw$rs6971 == 1 ] <- "AG"
raw$rs6971[raw$rs6971 == 2 ] <- "GG"
table(raw$rs6971)

df$rs6971_WES <- raw$rs6971[match(df$SJLID, raw$IID)]

write.table(df, "Z:/ResearchHome/ClusterHome/aneupane/data/request/Kyla/Within_50m_genotype_AN.txt", col.names =  T, sep = "\t", row.names = F)
