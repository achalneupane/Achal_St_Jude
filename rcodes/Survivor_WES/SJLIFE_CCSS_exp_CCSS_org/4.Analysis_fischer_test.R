library(data.table)
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife/phenotype.RData")
## 1. clinvar
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar <- clinvar[!duplicated(clinvar$SNP),]
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)

# MAF <- 0.01
MAF <- 0.0001

## Raw file maf
MAFname <- format(as.numeric(MAF), scientific = FALSE, digits = 5)
MAFname <- "" # no MAF filter in raw

# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar.all <- clinvar[which(clinvar$AF < MAF),]
clinvar.eur <- clinvar[which(clinvar$AF < MAF & clinvar$AF_nfe <  MAF),]
clinvar.afr <- clinvar[which(clinvar$AF < MAF & clinvar$AF_afr <  MAF),]


loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)
# - 
# 70494 

loftee <- loftee[!duplicated(loftee$SNP),]


loftee$AF <- as.numeric(loftee$AF.1)
loftee$AF_nfe <- as.numeric(loftee$AF_nfe)
loftee$AF_afr <- as.numeric(loftee$AF_afr)
loftee$CHROM  <- sub("([0-9XY]+):.+", "\\1", loftee$SNP)
loftee$POS <- sub("chr[0-9XY]+:(\\d+):.+", "\\1", loftee$SNP)

# make rare
loftee.all <- loftee[which(loftee$AF < MAF),]
loftee.eur <- loftee[which(loftee$AF < MAF & loftee$AF_nfe <  MAF),]
loftee.afr <- loftee[which(loftee$AF < MAF & loftee$AF_afr <  MAF),]


snpeff <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/missense_variants_with_overlap_in_more_than_90_percent_of_prediction_tools_all_cols.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(snpeff)
snpeff$SNP <- sub(";.*", "", snpeff$ID)
snpeff <- snpeff[!duplicated(snpeff$SNP),]

snpeff$AF <- as.numeric(snpeff$AF)
snpeff$AF_nfe <- as.numeric(snpeff$AF_nfe)
snpeff$AF_afr <- as.numeric(snpeff$AF_afr)


# make rare
snpeff.all <- snpeff[which(snpeff$AF < MAF),]
snpeff.eur <- snpeff[which(snpeff$AF < MAF & snpeff$AF_nfe <  MAF),]
snpeff.afr <- snpeff[which(snpeff$AF < MAF & snpeff$AF_afr <  MAF),]


cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
# 46076 
# 
cc <- c(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))

bim.QC.sjlife.PLP <- fread(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant//sjlife/all_rare_variants", MAFname,"_all_sjlife.bim"))
table(cc %in% bim.QC.sjlife.PLP$V2)
# FALSE  TRUE 
# 16818 29254 # maf 0.01
# 34256 11820 # maf 0.0001

raw <- fread(paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant//sjlife/all_rare_variants", MAFname,"_all_sjlife_recodeA.raw"))
sum(duplicated(raw$IID))
# 0
dim(raw)
# 4516 11826 # maf 0.0001
# 4516 29406 no maf filtered raw
# raw <- raw[!(duplicated(raw$IID)),]
raw <- as.data.frame(raw)
rownames(raw) <- raw$IID
raw <- raw[,-c(1:6)]
# HEADER <- colnames(raw)[-c(1:6)]
HEADER=colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(bim.QC.sjlife.PLP$V2 %in% colnames(raw))
# TRUE 
# 29254 # maf 0.01
# 11820 # maf 0.0001
# 29400 # no maf filtered raw
# colnames(raw)[!colnames(raw) %in% bim.QC.sjlife.PLP$V2]
raw$IID <- rownames(raw)

## run fisher exact test with grade 2 or higher ones with clinvar.all
raw.clinvar.all <- raw[colnames(raw) %in% c("IID", clinvar.all$SNP) ]
dim(raw.clinvar.all)
# 4516 6452 # gnomAD maf 0.01
# 4516 5276 # gnomAD maf 0.0001

raw.clinvar.all$carrier <- ifelse(rowSums(raw.clinvar.all[grepl("chr", colnames(raw.clinvar.all))], na.rm = T) > 0, 1, 0)
table(raw.clinvar.all$carrier)
# 0    1 
# 598 3918 
# 1698 2818 ## gnomAD maf 0.0001

CTCAE.data.4 <- CTCAE.data.4[CTCAE.data.4$sjlid %in% raw.clinvar.all$IID,]
CTCAE.data.4$carrier <- raw.clinvar.all$carrier[match(CTCAE.data.4$sjlid, raw.clinvar.all$IID)]

table(CTCAE.data.4$carrier)

##########################################################
## run fisher's exact test for grade greater than zero! ##
##########################################################
# Create an empty dataframe to store the results
results <- data.frame(variable = character(),
                      pvalue = numeric(),
                      OR.CI = character(),
                      stringsAsFactors = FALSE)

# Loop over all columns that end with '_status_gt_0'
for (col_name in names(CTCAE.data.4)) {
  if (grepl("_status_gt_0$", col_name)) {
    
    # Try to apply Fisher's test for each column with carrier
    tryCatch({
      gene.test <- fisher.test(table(CTCAE.data.4[[col_name]], CTCAE.data.4$carrier))
      pvalue <- gene.test$p.value
      OR.CI <- paste0(round(gene.test$estimate, 2), " (", 
                      paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
    }, error = function(e) {
      pvalue <- NA
      OR.CI <- NA
    })
    
    # Append the result to the dataframe
    results <- rbind(results, data.frame(variable = col_name,
                                         pvalue = pvalue,
                                         OR.CI = OR.CI,
                                         stringsAsFactors = FALSE))
  }
}

# View the final results
print(results)


#########################################################
## run fisher's exact test for grade greater than two! ##
#########################################################
# Create an empty dataframe to store the results
results <- data.frame(variable = character(),
                      pvalue = numeric(),
                      OR.CI = character(),
                      stringsAsFactors = FALSE)

# Loop over all columns that end with '_status_gt_2'
for (col_name in names(CTCAE.data.4)) {
  if (grepl("_status_gt_2$", col_name)) {
    
    # Try to apply Fisher's test for each column with carrier
    tryCatch({
      gene.test <- fisher.test(table(CTCAE.data.4[[col_name]], CTCAE.data.4$carrier))
      pvalue <- gene.test$p.value
      OR.CI <- paste0(round(gene.test$estimate, 2), " (", 
                      paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
    }, error = function(e) {
      pvalue <- NA
      OR.CI <- NA
    })
    
    # Append the result to the dataframe
    results <- rbind(results, data.frame(variable = col_name,
                                         pvalue = pvalue,
                                         OR.CI = OR.CI,
                                         stringsAsFactors = FALSE))
  }
}

# View the final results
print(results)


###########################################################
## run fisher's exact test for grade greater than three! ##
###########################################################
# Create an empty dataframe to store the results
results <- data.frame(variable = character(),
                      pvalue = numeric(),
                      OR.CI = character(),
                      stringsAsFactors = FALSE)

# Loop over all columns that end with '_status_gt_3'
for (col_name in names(CTCAE.data.4)) {
  if (grepl("_status_gt_3$", col_name)) {
    
    # Try to apply Fisher's test for each column with carrier
    tryCatch({
      gene.test <- fisher.test(table(CTCAE.data.4[[col_name]], CTCAE.data.4$carrier))
      pvalue <- gene.test$p.value
      OR.CI <- paste0(round(gene.test$estimate, 2), " (", 
                      paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
    }, error = function(e) {
      pvalue <- NA
      OR.CI <- NA
    })
    
    # Append the result to the dataframe
    results <- rbind(results, data.frame(variable = col_name,
                                         pvalue = pvalue,
                                         OR.CI = OR.CI,
                                         stringsAsFactors = FALSE))
  }
}

# View the final results
print(results)
