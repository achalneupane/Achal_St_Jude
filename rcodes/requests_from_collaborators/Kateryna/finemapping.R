## Chunk 1.
# CMP_chr1_22.assoc.logistic.clean.Psorted
gwas.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb", sep = " ", header = T)
# gwas.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Summary_results_AN.txt", sep = "\t", header = T)
dim(gwas.dat)
gwas.dat.original <- gwas.dat 
# gwas.dat$SNP <- paste(gwas.dat$CHR, gwas.dat$BP, gwas.dat$A2, gwas.dat$A1, sep = ":")
# gwas.dat$SNP <- gsub("chr", "", gwas.dat$SNP)
# sjlife.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife.bim", sep = "\t", header = F)

## Depends on which data you want to use. QCed or without QC
# sjlife.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_geno.0.1_hwe.1e-10.bim", sep = "\t", header = F) # QCed
sjlife.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed.bim", sep = "\t", header = F)
dim(sjlife.bim)
# sjlife.bim$KEY <- paste(sjlife.bim$V1, sjlife.bim$V4, sjlife.bim$V6, sjlife.bim$V5, sep = ":")
sum(gwas.dat$SNP %in% sjlife.bim$V2)
# 827

gwas.dat$SNP[(!gwas.dat$SNP %in% sjlife.bim$V2)]

# sjlife.bim$V2 <- gsub("chr", "", sjlife.bim$V2)

sjlife.bim <- sjlife.bim[sjlife.bim$V2 %in% gwas.dat$SNP,]
sjlife.bim$V2[duplicated(sjlife.bim$KEY)]

dim(sjlife.bim)
# [1] 859    6

# wanted.snp.list <- paste0("chr", sjlife.bim$V2)

# Merge sjlife bim with summary stat
gwas.dat <- cbind.data.frame(sjlife.bim, gwas.dat[match(sjlife.bim$V2, gwas.dat$SNP),])

sum(gwas.dat$V5 == gwas.dat$A1)
# 679
sum(gwas.dat$V6 == gwas.dat$A1)
# 180
# 679+180 = 859

# freq.dat <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_freq_993samples.frq_edited1", sep = "\t")
freq.dat <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.sjlid_chr16.PASS.decomposed_freq.frq_edited1", sep = "\t")
dim(freq.dat)

## Add MAF
sum(gwas.dat$V2 %in% freq.dat$SNP)
# 859 

gwas.dat$maf <- freq.dat$MAF[match(gwas.dat$V2, freq.dat$SNP)]
gwas.dat$maf.allele <- freq.dat$A1[match(gwas.dat$V2, freq.dat$SNP)]
gwas.dat$maf[gwas.dat$A1 != gwas.dat$maf.allele] <- 1-gwas.dat$maf[gwas.dat$A1 != gwas.dat$maf.allele]


gwas.dat$lowerCI <-  as.numeric(gwas.dat$L95)
gwas.dat$upperCI <-  as.numeric(gwas.dat$U95)


gwas.dat$beta <- log(gwas.dat$OR)
gwas.dat$lowerCI.beta <- log(gwas.dat$lowerCI)
gwas.dat$upperCI.beta <- log(gwas.dat$upperCI)

## upper C.I of beta = beta + se(beta) x 1.96
# gwas.dat$beta.SE <- (gwas.dat$upperCI.beta - gwas.dat$beta)/1.96
gwas.dat$beta.SE <- as.numeric(gwas.dat$SE) ## SE remains same
## Get A1 and A2
gwas.dat$allele1 <- gwas.dat$A1
gwas.dat$allele2 <- gwas.dat$A1

gwas.dat$allele2[gwas.dat$A1!=gwas.dat$V6] <- gwas.dat$V6[gwas.dat$A1!=gwas.dat$V6]
gwas.dat$allele2[gwas.dat$A1!=gwas.dat$V5] <- gwas.dat$V5[gwas.dat$A1!=gwas.dat$V5]

# Now, finally flip the allele and change the sign of beta where maf > 0.5 (50%)
gwas.dat[gwas.dat$maf > 0.5 , c("allele1", "allele2")] <- gwas.dat[gwas.dat$maf > 0.5, c("allele2", "allele1")] # swap alleles
gwas.dat$beta[gwas.dat$maf > 0.5] <- (-1 * gwas.dat$beta[gwas.dat$maf > 0.5]) # Change the sign of beta
gwas.dat$maf[gwas.dat$maf > 0.5] <- 1-gwas.dat$maf[gwas.dat$maf > 0.5] # Finally, change the maf values


## Now Keep only the the required columns
gwas.dat <- cbind.data.frame(rsid=gwas.dat$SNP, chromosome=gwas.dat$V1, position=gwas.dat$V4, 
                             allele1 = gwas.dat$allele1, allele2 = gwas.dat$allele2, maf = gwas.dat$maf, beta = gwas.dat$beta, se = gwas.dat$beta.SE, P=gwas.dat$P)


# # Write data.z file for FINEMAP
sum(gwas.dat$maf < 0.01)
# 0

## First only keep those with maf >1 %; then separate the summary statistics for chr16 
gwas.dat.maf.gt.1perc <- gwas.dat[!gwas.dat$maf < 0.01,]

gwas.dat.maf.gt.1perc <- gwas.dat.maf.gt.1perc[grepl("16", gwas.dat.maf.gt.1perc$chromosome),]
write.table(as.data.frame(gwas.dat.maf.gt.1perc$rsid), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap_kateryna/samplesnp_gt_MAF_1_perc_vars_meta.list", quote = F, row.names = F, col.names = F, sep = " ")
# After extracting these variants from plink, I am re-ordering the summary stat
# to match the order of varinats in plink. This is to ensure LD and summary stat
# have same variant order


samplesnp_gt_MAF_1_perc_vars.dat.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap_kateryna/chr16_finemap_plink.bim")
gwas.dat.maf.gt.1perc <- gwas.dat.maf.gt.1perc[match(samplesnp_gt_MAF_1_perc_vars.dat.bim$V2, gwas.dat.maf.gt.1perc$rsid),]
table(gwas.dat.maf.gt.1perc$rsid == samplesnp_gt_MAF_1_perc_vars.dat.bim$V2)
# TRUE 
# 859 

write.table(gwas.dat.maf.gt.1perc, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap_kateryna/samplesnp_gt_MAF_1_perc.z", quote = F, row.names = F, col.names = T, sep = " ")


## Input summary stat file for COJO analysis 
chr16.COJO <- cbind.data.frame(SNP=gwas.dat.maf.gt.1perc$rsid, A1=gwas.dat.maf.gt.1perc$allele1, A2=gwas.dat.maf.gt.1perc$allele2, 
                              freq=gwas.dat.maf.gt.1perc$maf, b=gwas.dat.maf.gt.1perc$beta, se=gwas.dat.maf.gt.1perc$se, p=gwas.dat.maf.gt.1perc$P)
chr16.COJO$N <- 993

## COJO files
write.table(chr16.COJO, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap_kateryna/samplesnp_gt_MAF_1_perc_vars.ma", quote = F, row.names = F, col.names = F, sep = " ")
