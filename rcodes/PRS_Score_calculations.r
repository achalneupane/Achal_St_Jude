## Michigan overall; Take POS_START
OVERALL.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/OVERALL_PRSWEB_PHECODE174.1_Onco-iCOGS-Overall-BRCA_PRS-CS_UKB_20200608_WEIGHTS_edited_1.txt", header = T)
OVERALL.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_OVERALL_GrCh38.bed", header = F)
colnames(OVERALL.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

OVERALL.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(OVERALL.38[match(OVERALL.37$`ChrCHROM:POS`, OVERALL.38$KEY),], OVERALL.37)

FINAL.dat <- cbind(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
write.table(FINAL.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_OVERALL_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")



## Michigan ER_NEG; Take POS_START
ER_NEG.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/ER_NEG_PRSWEB_PHECODE174.1_GWAS-Catalog-r2019-05-03-X174.1_PT_UKB_20200608_WEIGHTS_edited_1.txt", header = T)
ER_NEG.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_NEG_GrCh38.bed", header = F)
colnames(ER_NEG.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

ER_NEG.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(ER_NEG.38[match(ER_NEG.37$`ChrCHROM:POS`, ER_NEG.38$KEY),], ER_NEG.37)

FINAL.dat <- cbind(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
write.table(FINAL.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_NEG_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")


## Michigan ER_POS; Take POS_START
ER_POS.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/ER_POS_PRSWEB_PHECODE174.1_Onco-iCOGS-ER-positive-BRCA_PRS-CS_MGI_20200608_WEIGHTS_edited_1.txt", header = T)
ER_POS.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_POS_GrCh38.bed", header = F)
colnames(ER_POS.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

ER_POS.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(ER_POS.38[match(ER_POS.37$`ChrCHROM:POS`, ER_POS.38$KEY),], ER_POS.37)

FINAL.dat <- cbind(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
write.table(FINAL.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_POS_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")


