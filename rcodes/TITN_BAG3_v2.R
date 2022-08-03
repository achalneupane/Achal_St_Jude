############################
## TITN and BAG3 analysis ##
############################

## Chunk 1.
gwas.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Summary_results_AN.txt", sep = "\t", header = T)
dim(gwas.dat)
gwas.dat$SNP <- paste(gwas.dat$SNP, gwas.dat$A2, gwas.dat$A1, sep = ":")


sjlife.bim <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife.bim", sep = "\t", header = F)
dim(sjlife.bim)
sjlife.bim$KEY <- paste(sjlife.bim$V1, sjlife.bim$V4, sjlife.bim$V6, sjlife.bim$V5, sep = ":")


sum(gwas.dat$SNP %in% sjlife.bim$KEY)
# 2601

sjlife.bim$V2 <- gsub("chr", "", sjlife.bim$V2)

sjlife.bim <- sjlife.bim[sjlife.bim$KEY %in% gwas.dat$SNP,]
sjlife.bim$V2[duplicated(sjlife.bim$KEY)]

## checked dbSNP and these seem to be the real ones
# 10:119588327:TTTTC:T
# 10:119427920:ATCACTGCCATCATCAGCC:A # matches A1 and A2 as in summary stat
# 10:119485655:A:ATATTTTATTTTATTT
# 10:119528844:T:TAATA; maf 0.07029 # matches A1 and A2 as in summary stat
# 10:119614280:A:AAAAT; 
# 10:119614280:A:AAAAT;10:119614280:A:AAAATAAATAAAT     AAAAT     A       0.1388  3270
## There are variants 



# Based on the conversation with Yadav on 08/03/2022, I am deleting 10:119614280
# and 10:119485655 from the list. I am also going to swap alleles for
# 10:119427920 so that minor alleles correspond to the effect estimate.

# In conclusion,
# 10:119427920:ATCACTGCCATCATCAGCC:A # Keep and swap alleles
# 10:119427920:A:ATCACTGCCATCATCAGCC # remove
# 10:119528844:T:TAATA # keep
# 10:119528844:TAATA:T # remove
# 10:119588327:TTTTC:T # keep
# 10:119588327:T:TTTTC # remove
# 10:119614280 # remove
# 10:119485655 # remove



sjlife.bim$V2[grepl("10:119528844|10:119427920|10:119588327", sjlife.bim$V2)]
# [1] "10:119427920:A:ATCACTGCCATCATCAGCC" "10:119427920:ATCACTGCCATCATCAGCC:A" "10:119528844:T:TAATA"               "10:119528844:TAATA:T"              
# [5] "10:119588327:T:TTTTC"               "10:119588327:TTTTC:T"  

sjlife.bim <- sjlife.bim[!grepl("10:119614280|10:119485655|10:119588327:T:TTTTC|10:119528844:TAATA:T|10:119427920:A:ATCACTGCCATCATCAGCC", sjlife.bim$V2),]
# dim(sjlife.bim)
# [1] 2599    7

wanted.snp.list <- paste0("chr", sjlife.bim$V2)

write.table(as.data.frame(wanted.snp.list), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.list", quote = F, row.names = F, col.names = F)

# Merge sjlife bim with summary stat
gwas.dat <- cbind.data.frame(sjlife.bim, gwas.dat[match(sjlife.bim$KEY, gwas.dat$SNP),])

sum(gwas.dat$V5 == gwas.dat$A1)
# 2599
sum(gwas.dat$V6 == gwas.dat$A2)
# 2599

freq.dat <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/sjlife.freq.out.frq_edited1", sep = "\t")
dim(freq.dat)
freq.dat


# freq.dat$KEY <- vapply(strsplit(freq.dat$SNP, ':'), function(x)  paste(x[seq.int(2)], collapse=':'), character(1L))
# freq.dat$KEY <- paste(freq.dat$KEY, freq.dat$A2, freq.dat$A1, sep = ":")


## Add MAF
sum(gwas.dat$V2 %in% freq.dat$SNP)
# 2599

gwas.dat$maf <- freq.dat$MAF[match(gwas.dat$V2, freq.dat$SNP)]



gwas.dat <- cbind.data.frame(rsid=gwas.dat$V2, SNP=gwas.dat$SNP, chromosome=gwas.dat$CHR, position=gwas.dat$BP, 
                             allele1 = gwas.dat$A1, allele2 = gwas.dat$A2, maf = gwas.dat$maf, OR_95percent_CI_SJLIFE = gwas.dat$OR_95percent_CI_SJLIFE)

gwas.dat$OR <- as.numeric(sapply(strsplit(gwas.dat$OR_95percent_CI_SJLIFE, " "), `[`, 1))

gwas.dat$CI <- gsub("\\(|\\)| ", "",sapply(strsplit(gwas.dat$OR_95percent_CI_SJLIFE, " "), `[`, 2))
gwas.dat$lowerCI <-  as.numeric(sapply(strsplit(gwas.dat$CI, "-"), `[`, 1))
gwas.dat$upperCI <-  as.numeric(sapply(strsplit(gwas.dat$CI, "-"), `[`, 2))


gwas.dat$beta <- log(gwas.dat$OR)
gwas.dat$lowerCI.beta <- log(gwas.dat$lowerCI)
gwas.dat$upperCI.beta <- log(gwas.dat$upperCI)

## upper C.I of beta = beta + se(beta) x 1.96
gwas.dat$beta.SE <- (gwas.dat$upperCI.beta - gwas.dat$beta)/1.96



# swap allele for this variant and change the sign of beta 
gwas.dat$rsid == "10:119427920:ATCACTGCCATCATCAGCC:A"

gwas.dat[gwas.dat$rsid == "10:119427920:ATCACTGCCATCATCAGCC:A", c("allele1", "allele2")] <- gwas.dat[gwas.dat$rsid == "10:119427920:ATCACTGCCATCATCAGCC:A", c("allele2", "allele1")]

# Also change the sign of beta
gwas.dat$beta[gwas.dat$rsid == "10:119427920:ATCACTGCCATCATCAGCC:A"]  <- (-1 * gwas.dat$beta[gwas.dat$rsid == "10:119427920:ATCACTGCCATCATCAGCC:A"])


# Now, finally flip the allele and change the sign of beta where maf > 0.5
gwas.dat[gwas.dat$maf > 0.5 , c("allele1", "allele2")] <- gwas.dat[gwas.dat$maf > 0.5, c("allele2", "allele1")] # swap alleles
gwas.dat$beta[gwas.dat$maf > 0.5] <- (-1 * gwas.dat$beta[gwas.dat$maf > 0.5]) # Change the sign of beta
gwas.dat$maf[gwas.dat$maf > 0.5] <- 1-gwas.dat$maf[gwas.dat$maf > 0.5] # Finally, change the maf values


## Now Keep only the the required columns
gwas.dat <- cbind.data.frame(rsid=gwas.dat$rsid, chromosome=gwas.dat$chromosome, position=gwas.dat$position, 
                             allele1 = gwas.dat$allele1, allele2 = gwas.dat$allele2, maf = gwas.dat$maf, beta = gwas.dat$beta, se = gwas.dat$beta.SE)


# Now sort this summary stat in the same order as Plink
gwas.dat$rsid <- paste0("chr", gwas.dat$rsid)
samplesnp.dat.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.dat.bim", sep = "\t", header = F)

sum(samplesnp.dat.bim$V2 %in% gwas.dat$rsid)
# 2599
gwas.dat <- gwas.dat[match(samplesnp.dat.bim$V2, gwas.dat$rsid),]
table(gwas.dat$rsid == samplesnp.dat.bim$V2)
# TRUE 
# 2599 

# Write data.z file for FINEMAP
write.table(gwas.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.z", quote = F, row.names = F, col.names = T, sep = " ")


