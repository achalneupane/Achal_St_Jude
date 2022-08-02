############################
## TITN and BAG3 analysis ##
############################

## Chunk 1.
gwas.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Summary_results_AN.txt", sep = "\t", header = T)
dim(gwas.dat)
freq.dat <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife.freq.out.frq_edited1", sep = "\t")
dim(freq.dat)

gwas.dat$KEY1 <- paste0(gwas.dat$SNP, ":", gwas.dat$A2, ":", gwas.dat$A1)
gwas.dat$KEY2 <- paste0(gwas.dat$SNP, ":", gwas.dat$A1, ":", gwas.dat$A2)

sum(gwas.dat$KEY1 %in% freq.dat$SNP)
# 2302
sum(gwas.dat$KEY2 %in% freq.dat$SNP)
# 316

gwas.dat$maf_same <- freq.dat$MAF[match(gwas.dat$KEY1, freq.dat$SNP)]
gwas.dat$maf_swapped <- freq.dat$MAF[match(gwas.dat$KEY2, freq.dat$SNP)]



gwas.dat$maf <- gwas.dat$maf_same

## SNP ID should be same as plink, so keeping the order of alleles accordingly
gwas.dat$rsid <- gwas.dat$KEY1
gwas.dat$rsid[is.na(gwas.dat$maf)] <- gwas.dat$KEY2[is.na(gwas.dat$maf)]

gwas.dat$allele_swapped <- ifelse(is.na(gwas.dat$maf), "Y", "N")



# Now adding maf of SNPs with swapped alleles
gwas.dat$maf[gwas.dat$allele_swapped == "Y"] <- 1-gwas.dat$maf_swapped[gwas.dat$allele_swapped == "Y"]
sum(is.na(gwas.dat$maf))
# 0

gwas.dat <- cbind.data.frame(rsid=gwas.dat$rsid, chromosome=gwas.dat$CHR, position=gwas.dat$BP, 
                             allele1 = gwas.dat$A1, allele2 = gwas.dat$A2, maf = gwas.dat$maf, OR_95percent_CI_SJLIFE = gwas.dat$OR_95percent_CI_SJLIFE, allele_swapped = gwas.dat$allele_swapped)

gwas.dat$OR <- as.numeric(sapply(strsplit(gwas.dat$OR_95percent_CI_SJLIFE, " "), `[`, 1))

gwas.dat$CI <- gsub("\\(|\\)| ", "",sapply(strsplit(gwas.dat$OR_95percent_CI_SJLIFE, " "), `[`, 2))
gwas.dat$lowerCI <-  as.numeric(sapply(strsplit(gwas.dat$CI, "-"), `[`, 1))
gwas.dat$upperCI <-  as.numeric(sapply(strsplit(gwas.dat$CI, "-"), `[`, 2))


gwas.dat$beta <- log(gwas.dat$OR)
gwas.dat$lowerCI.beta <- log(gwas.dat$lowerCI)
gwas.dat$upperCI.beta <- log(gwas.dat$upperCI)

## upper C.I of beta = beta + se(beta) x 1.96
gwas.dat$beta.SE <- (gwas.dat$upperCI.beta - gwas.dat$beta)/1.96


## Now Keep only the the required columns
gwas.dat <- cbind.data.frame(rsid=gwas.dat$rsid, chromosome=gwas.dat$chromosome, position=gwas.dat$position, 
                             allele1 = gwas.dat$allele1, allele2 = gwas.dat$allele2, maf = gwas.dat$maf, beta = gwas.dat$beta, se = gwas.dat$beta.SE)


gwas.dat$rsid <- paste0("chr",gwas.dat$rsid)


## write SNP list to get LD-matrix
write.table(as.data.frame(gwas.dat$rsid), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.list", quote = F, row.names = F, col.names = F)


# Error : Expected a floating point value greater than 0.0 and smaller or equal than 0.5 in line 3 column 6 of file 'samplesnp.z'!
gwas.dat.edited <- gwas.dat  
gwas.dat.edited$maf

write.table(gwas.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.z", quote = F, row.names = F, col.names = T, sep = " ")


## Reorder LD file
ld.file <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp_ld.ld", header = F)
dim(ld.file)

sum(ld.file$SNP_B %in% gwas.dat$rsid)
