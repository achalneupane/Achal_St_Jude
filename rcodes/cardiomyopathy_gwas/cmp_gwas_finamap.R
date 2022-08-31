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



freq.dat <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap/sjlife.freq.out.frq_edited1", sep = "\t")
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
# gwas.dat[gwas.dat$maf > 0.5 , c("allele1", "allele2")] <- gwas.dat[gwas.dat$maf > 0.5, c("allele2", "allele1")] # swap alleles
# gwas.dat$beta[gwas.dat$maf > 0.5] <- (-1 * gwas.dat$beta[gwas.dat$maf > 0.5]) # Change the sign of beta
# gwas.dat$maf[gwas.dat$maf > 0.5] <- 1-gwas.dat$maf[gwas.dat$maf > 0.5] # Finally, change the maf values


# Now sort this summary stat in the same order as Plink
samplesnp.dat.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap/samplesnp.dat.bim", sep = "\t", header = F)

sum(samplesnp.dat.bim$V2 %in% gwas.dat$rsid)
# 859
gwas.dat <- gwas.dat[match(samplesnp.dat.bim$V2, gwas.dat$rsid),]
table(gwas.dat$rsid == samplesnp.dat.bim$V2)
# TRUE 
# 859 


# The GWAS was run on maf 0.05
sum(gwas.dat$maf < 0.05)
# 0

# Write data.z file for FINEMAP
write.table(gwas.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap/samplesnp.z", quote = F, row.names = F, col.names = T, sep = " ")











## Input summary stat file for COJO analysis 
TITN.COJO <- cbind.data.frame(SNP=gwas.dat.maf.gt.1perc.TITN$rsid, A1=gwas.dat.maf.gt.1perc.TITN$allele1, A2=gwas.dat.maf.gt.1perc.TITN$allele2, 
                              freq=gwas.dat.maf.gt.1perc.TITN$maf, b=gwas.dat.maf.gt.1perc.TITN$beta, se=gwas.dat.maf.gt.1perc.TITN$se, p=gwas.dat.maf.gt.1perc.TITN$P)
TITN.COJO$N <- 1645

BAG3.COJO <- cbind.data.frame(SNP=gwas.dat.maf.gt.1perc.BAG3$rsid, A1=gwas.dat.maf.gt.1perc.BAG3$allele1, A2=gwas.dat.maf.gt.1perc.BAG3$allele2, 
                              freq=gwas.dat.maf.gt.1perc.BAG3$maf, b=gwas.dat.maf.gt.1perc.BAG3$beta, se=gwas.dat.maf.gt.1perc.BAG3$se, p=gwas.dat.maf.gt.1perc.BAG3$P)
BAG3.COJO$N <- 1645

## COJO files
write.table(TITN.COJO, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test/samplesnp_TITN_gt_MAF_1_perc_vars.ma", quote = F, row.names = F, col.names = F, sep = " ")
write.table(BAG3.COJO, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/cojo_test/samplesnp_BAG3_gt_MAF_1_perc_vats.ma", quote = F, row.names = F, col.names = F, sep = " ")



#########################
## Reading the results ##
#########################
# Assuming 1 causal
samplesnp_TITN_gt_MAF_1_perc.snp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP_results/assuming_1_causal_SNP/samplesnp_TITN_gt_MAF_1_perc.snp", header = T)
dim(samplesnp_TITN_gt_MAF_1_perc.snp)
# 947
max(samplesnp_TITN_gt_MAF_1_perc.snp$prob)
# 0.014; 2:178562809 is top 27th

samplesnp_BAG3_gt_MAF_1_perc.snp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP_results/assuming_1_causal_SNP/samplesnp_BAG3_gt_MAF_1_perc.snp", header = T)
dim(samplesnp_BAG3_gt_MAF_1_perc.snp)
max(samplesnp_BAG3_gt_MAF_1_perc.snp$prob)
# 0.038; 10:119670121 is top 74th


prob.stats.titn.bag3.maf.gt.1.perc <- rbind.data.frame(samplesnp_TITN_gt_MAF_1_perc.snp, samplesnp_BAG3_gt_MAF_1_perc.snp)
dim(prob.stats.titn.bag3.maf.gt.1.perc)

prob.stats.titn.bag3.maf.gt.1.perc$SNP <- paste(prob.stats.titn.bag3.maf.gt.1.perc$chromosome, prob.stats.titn.bag3.maf.gt.1.perc$position, sep = ":")

sum(gwas.dat.original$SNP %in% prob.stats.titn.bag3.maf.gt.1.perc$SNP)
# 2541

gwas.dat.original$FINEMAP_Prob <- prob.stats.titn.bag3.maf.gt.1.perc$prob[match(gwas.dat.original$SNP, prob.stats.titn.bag3.maf.gt.1.perc$SNP)]
gwas.dat$SNP <- paste(gwas.dat$chromosome, gwas.dat$position, sep = ":")
gwas.dat.original$MAF <- gwas.dat$maf[match(gwas.dat.original$SNP, gwas.dat$SNP)]

write.table(gwas.dat.original, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP_prob.stats.titn.bag3.maf.gt.1.perc.txt", quote = F, row.names = F, col.names = T, sep = "\t")


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/FINEMAP.RDATA")


