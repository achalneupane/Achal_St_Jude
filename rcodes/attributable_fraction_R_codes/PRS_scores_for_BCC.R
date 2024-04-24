## Prepare files
# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia/PRS_CCSS/")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//attr_fraction/prs/NMSC")
library(data.table)
library(dplyr)
#######################
## PGS000688 _GRCh38 ##
#######################

PGS000119.grch37 <- read.delim("PGS000119_GRCh37", sep = "\t", header = T)
PGS000119.grch38 <- read.delim("PGS000119", sep = "\t", header = T)

# get overlapping variants
ccss.org.bim <- read.table("CCSS_org_bim_vars_PGS000119")

replace.allele <- c("G", "A", "C", "G", "T", "C", "A", "C", "A", "A", "A", "T", "A", "G", "G", "T", "T", "T", "C", "A", "T", "T", "G", "C", "C", "C", "C", "T", "A", "T", "G", "T")
 
PGS000119.grch37$other_allele <- replace.allele
PGS000119.grch38$other_allele <- replace.allele

PGS000119.grch37 <- PGS000119.grch37[,c("hm_chr", "hm_pos", "other_allele", "effect_allele", "effect_weight")]
colnames(PGS000119.grch37) <- c("CHROM", "POS_GRCh37", "REF", "Effect_allele", "Effect_size")

PGS000119.grch38 <- PGS000119.grch38[,c("hm_chr", "hm_pos", "other_allele", "effect_allele", "effect_weight")]
colnames(PGS000119.grch38) <- c("CHROM", "POS_GRCh38", "REF", "Effect_allele", "Effect_size")

################################
BCC <- PGS000119.grch38

## Swap alleles if beta is negative
BCC[BCC$Effect_size < 0, c("REF", "Effect_allele")] <- BCC[BCC$Effect_size < 0, c("Effect_allele", "REF")]
BCC$Effect_size <- abs(BCC$Effect_size)

BCC$POS_GRCh38 <- as.numeric(BCC$POS_GRCh38)
BCC$CHROM <- paste0("chr", BCC$CHROM)

write.table(BCC, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//attr_fraction/prs/NMSC/PGS000119.grch38.dat", col.names = F, row.names = F, quote = F)

## GRCh37
BCC <- PGS000119.grch37

## Swap alleles if beta is negative
BCC[BCC$Effect_size < 0, c("REF", "Effect_allele")] <- BCC[BCC$Effect_size < 0, c("Effect_allele", "REF")]
BCC$Effect_size <- abs(BCC$Effect_size)

BCC$POS_GRCh37 <- as.numeric(BCC$POS_GRCh37)
BCC$CHROM <- paste0("chr", BCC$CHROM)

write.table(BCC, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//attr_fraction/prs/NMSC/PGS000119.grch37.dat", col.names = F, row.names = F, quote = F)

