############################
## TITN and BAG3 analysis ##
############################

## Chunk 1.
gwas.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb", sep = " ", header = T)
dim(gwas.dat)
head(gwas.dat)
gwas.dat.original <- gwas.dat 

sjlife.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap/sjlife_CMP.bim", sep = "\t", header = F)
dim(sjlife.bim)

gwas.dat$SNP %in% sjlife.bim$V2
sum(gwas.dat$SNP %in% sjlife.bim$V2)
# 859 # all match


PHENO=read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/sjlife_eur_dox_only_pcs.pheno", header = T)



# freq.dat <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap/sjlife.freq.out.frq_edited1", sep = "\t")
freq.dat <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap/samplesnp.dat.freq.out.frq_edited1", sep = "\t")
dim(freq.dat)
freq.dat


## Add MAF
sum(gwas.dat$SNP %in% freq.dat$SNP)
# 859

gwas.dat$maf <- freq.dat$MAF[match(gwas.dat$SNP, freq.dat$SNP)]

gwas.dat$freqA1 <- freq.dat$A1[match(gwas.dat$SNP, freq.dat$SNP)]
gwas.dat$freqA2 <- freq.dat$A2[match(gwas.dat$SNP, freq.dat$SNP)]
table(gwas.dat$freqA1 == gwas.dat$A1)
# 859 # all true; so gwas.dat$freqA2 is gwas.dat$A2
gwas.dat$A2 <- gwas.dat$freqA2


#######################################
## calculate standard error for beta ##
#######################################
## upper C.I of beta = beta + se(beta) x 1.96
gwas.dat$lowerCI.OR <- gwas.dat$OR - (gwas.dat$SE)*1.96
gwas.dat$upperCI.OR <- gwas.dat$OR + (gwas.dat$SE)*1.96

gwas.dat$lowerCI.beta <- log(gwas.dat$lowerCI.OR)
gwas.dat$upperCI.beta <- log(gwas.dat$upperCI.OR)

gwas.dat$beta <- log(gwas.dat$OR)
## upper C.I of beta = beta + se(beta) x 1.96
gwas.dat$beta.SE <- (gwas.dat$upperCI.beta - gwas.dat$beta)/1.96


gwas.dat <- cbind.data.frame(rsid=gwas.dat$SNP, chromosome=gwas.dat$CHR, position=gwas.dat$BP, 
                             allele1 = gwas.dat$A1, allele2 = gwas.dat$A2, maf = gwas.dat$maf, beta = gwas.dat$beta, se=gwas.dat$beta.SE)


# # Now, finally flip the allele and change the sign of beta where maf > 0.5
gwas.dat[gwas.dat$maf > 0.5 , c("allele1", "allele2")] <- gwas.dat[gwas.dat$maf > 0.5, c("allele2", "allele1")] # swap alleles
gwas.dat$beta[gwas.dat$maf > 0.5] <- (-1 * gwas.dat$beta[gwas.dat$maf > 0.5]) # Change the sign of beta
gwas.dat$maf[gwas.dat$maf > 0.5] <- 1-gwas.dat$maf[gwas.dat$maf > 0.5] # Finally, change the maf values



# Filter for maf 0.01
sum(gwas.dat$maf < 0.01)
# 7

gwas.dat <- gwas.dat[!gwas.dat$maf < 0.01,]


# Yadav: I think there are some SNPs which are probably <1% in controls but higher in
# cases and this the overall freq is higher than 1%. Can you remove variants
# from the summary statistics file that have <1% in controls only and then
# condition for chr16:25611595:C:T only?

table(PHENO$CMP2plus)
# 1   2 
# 783 210 

CASES <- PHENO[PHENO$CMP2plus == 2, c("FID","IID")]
CONTROLS <- PHENO[PHENO$CMP2plus == 1, c("FID","IID")]

write.table(CONTROLS, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2/controls.list", quote = F, row.names = F, col.names = T, sep = " ")

# remove variants that are less than 1% in controls
controls.freq <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2/samplesnp.dat.controls.freq.out.frq_edited1", header = T)
controls.freq.lt.1maf <- controls.freq[controls.freq$MAF < 0.01,]
sum(gwas.dat$rsid %in%  controls.freq.lt.1maf$SNP)
# 0 # these were 7 variants already removed when filtered for maf 0.01

# extract from plink
write.table(gwas.dat$rsid, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2/final.snp.list", row.names = F, quote = F)

# Now sort this summary stat in the same order as Plink
samplesnp.dat.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2/samplesnp.dat.final.bim", sep = "\t", header = F)

sum(samplesnp.dat.bim$V2 %in% gwas.dat$rsid)
# 859
gwas.dat <- gwas.dat[match(samplesnp.dat.bim$V2, gwas.dat$rsid),]
table(gwas.dat$rsid == samplesnp.dat.bim$V2)
# TRUE 
# 852 

# Write data.z file for FINEMAP
write.table(gwas.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2/samplesnp.z", quote = F, row.names = F, col.names = T, sep = " ")


## FINEMAP results are in excel spreadsheet: Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2/samplesnp_finemap_results.xlsx




## Input summary stat file for COJO analysis 
COJO <- cbind.data.frame(SNP=gwas.dat$rsid, A1=gwas.dat$allele1, A2=gwas.dat$allele2, 
                            freq=gwas.dat$maf, b=gwas.dat$beta, se=gwas.dat$se)

COJO$P <- gwas.dat.original$P[match(COJO$SNP, gwas.dat.original$SNP)]

COJO$N <- 852


## COJO files
write.table(COJO, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap2/samplesnp_vars.ma", quote = F, row.names = F, col.names = F, sep = " ")




