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

sjlife.bim <- sjlife.bim[sjlife.bim$KEY %in% gwas.dat$SNP,]
sjlife.bim$KEY[duplicated(sjlife.bim$KEY)]

write.table(as.data.frame(wanted.snp.list), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.list", quote = F, row.names = F, col.names = F)




freq.dat <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife.freq.out.frq_edited1", sep = "\t")
dim(freq.dat)
freq.dat




freq.dat$KEY <- vapply(strsplit(freq.dat$SNP, ':'), function(x)  paste(x[seq.int(2)], collapse=':'), character(1L))
freq.dat$KEY <- paste(freq.dat$KEY, freq.dat$A2, freq.dat$A1, sep = ":")



sum(gwas.dat$SNP %in% freq.dat$KEY)
# 2601

gwas.dat$maf_original <- freq.dat$MAF[match(gwas.dat$SNP, freq.dat$KEY)]


# ## SNP ID should be same as plink, so keeping the order of alleles accordingly
# gwas.dat$rsid <- gwas.dat$KEY1
# gwas.dat$rsid[is.na(gwas.dat$maf)] <- gwas.dat$KEY2[is.na(gwas.dat$maf)]
# 
# gwas.dat$allele_swapped <- ifelse(is.na(gwas.dat$maf), "Y", "N")



# # Now adding maf of SNPs with swapped alleles
# gwas.dat$maf[gwas.dat$allele_swapped == "Y"] <- 1-gwas.dat$maf_swapped[gwas.dat$allele_swapped == "Y"]
# sum(is.na(gwas.dat$maf))
# # 0

gwas.dat <- cbind.data.frame(rsid=gwas.dat$SNP, chromosome=gwas.dat$CHR, position=gwas.dat$BP, 
                             allele1 = gwas.dat$A1, allele2 = gwas.dat$A2, maf = gwas.dat$maf, OR_95percent_CI_SJLIFE = gwas.dat$OR_95percent_CI_SJLIFE, maf_original = gwas.dat$maf_original)

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
# write.table(as.data.frame(paste(gwas.dat$rsid, gwas.dat$allele2, gwas.dat$allele1, sep = ":")), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.list", quote = F, row.names = F, col.names = F)
write.table(as.data.frame(gwas.dat$rsid), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.list", quote = F, row.names = F, col.names = F)


# Error : Expected a floating point value greater than 0.0 and smaller or equal than 0.5 in line 3 column 6 of file 'samplesnp.z'!
gwas.dat.edited <- gwas.dat  
gwas.dat.edited$maf

write.table(gwas.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp.z", quote = F, row.names = F, col.names = T, sep = " ")


## Reorder LD file
ld.file <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/samplesnp_ld.ld", header = F)
dim(ld.file)

sum(ld.file$SNP_B %in% gwas.dat$rsid)
