library(data.table)
## Missense
all.missense <- as.data.frame(fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/missense_variants_with_overlap_in_more_than_90_percent_of_prediction_tools_all_cols.txt", header = T, stringsAsFactors = F))
## Loftee
loftee <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/loftee/loftee_lines_HC_all_cols.txt", header = T, stringsAsFactors = F)
loftee2 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, stringsAsFactors = F)
dim(loftee)
loftee <- as.data.frame(loftee[loftee$SNP %in% loftee2$`#Uploaded_variation`,])


# Main folder: /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes
clinvar <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/all_new_clinvar_P_LP.txt", header = T, stringsAsFactors = F)
clinvar$SNP <- sub(";.*", "", clinvar$ID)

all.missense$from <- "missense"
loftee$from <- "loftee"
clinvar$from <- "clinvar"
clinvar <- rbind.data.frame(all.missense, loftee, clinvar)

clinvar$ID <- sub(";.*", "", clinvar$ID)
clinvar <- clinvar[!duplicated(clinvar$ID),]
# cc <- cbind.data.frame(clinvar$ID, clinvar$AF, clinvar$AF_nfe, clinvar$AF_afr, clinvar$CLNSIG)

write.table(as.data.frame(clinvar$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity//common/gwas/rare_variants/geneBased_June_24_2024/extract_PLP.txt", col.names = F, quote = F, sep = "\t")
dim(clinvar)
# 26328
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)

clinvar <- clinvar[,!duplicated(colnames(clinvar))]

colnames(clinvar) <- gsub("ANN\\[\\*\\]\\.", "", colnames(clinvar))
#############
## EUR PLP ##
#############
clinvar.nfe <- clinvar[which(clinvar$AF < 0.05 & clinvar$AF_nfe < 0.05),]
dim(clinvar.nfe)
# 18645
## Read MAC data with 0.05 maf
EUR_MAC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/tmp_EUR_merged_data_plink_MAC_EUR.acount", header = T, stringsAsFactors = F, skip = 'CHROM')
## Selecting variants that have MAC of at least 3
EUR_MAC <- EUR_MAC[EUR_MAC$ALT_CTS >= 3,]
clinvar.nfe <- clinvar.nfe[clinvar.nfe$ID %in% EUR_MAC$ID,] # These are also MAF < 0.05 in our EUR data
dim(clinvar.nfe)
# 572



library(dplyr)
# Keeping only those with 2 or more variants in each gene 
# Count the occurrences of each gene
gene_counts <- clinvar.nfe %>%
  group_by(GENE) %>%
  summarise(count = n())

# Merge the counts back into the original data frame
clinvar.nfe <- clinvar.nfe %>%
  left_join(gene_counts, by = "GENE") %>%
  rename(gene.counts.NFE = count)
# cc <- cbind.data.frame(clinvar.nfe$ID, clinvar.nfe$AF, clinvar.nfe$AF_nfe, clinvar.nfe$AF_afr, clinvar.nfe$CLNSIG, clinvar.nfe$gene.counts.NFE, clinvar.nfe$ANN....GENE)
dim(clinvar.nfe)
# 572

clinvar.nfe <- clinvar.nfe[clinvar.nfe$gene.counts.NFE >= 2,]
dim(clinvar.nfe)
# 102

table(clinvar.nfe$from)

write.table(as.data.frame(clinvar.nfe$ID), "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/vars.clinvar.missense.loftee.eur.maf0.05.mac3.txt", col.names = F, row.names = F, quote = F)

## Write genesets
EUR <- cbind.data.frame(CHROM=clinvar.nfe$CHROM, GENE=clinvar.nfe$GENE, SNPID=clinvar.nfe$ID)
write.table(EUR, "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/extract-variants-final-EUR.GENE-SNP", sep = "\t", row.names = F, quote = F, col.names = F)



#############
## AFR PLP ##
#############
clinvar.afr <- clinvar[which(clinvar$AF < 0.05 & clinvar$AF_afr < 0.05),]
dim(clinvar.afr)
# 18582
## Read MAC data with 0.05 maf
AFR_MAC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/tmp_AFR_merged_data_plink_MAC_AFR.acount", header = T, stringsAsFactors = F, skip = 'CHROM')
## Selecting variants that have MAC of at least 3
AFR_MAC <- AFR_MAC[AFR_MAC$ALT_CTS >= 3,]
clinvar.afr <- clinvar.afr[clinvar.afr$ID %in% AFR_MAC$ID,] # These are also MAF < 0.05 in our AFR data
dim(clinvar.afr)
# 223

library(dplyr)
# Keeping only those with 2 or more variants in each gene 
# Count the occurrences of each gene
gene_counts <- clinvar.afr %>%
  group_by(GENE) %>%
  summarise(count = n())

# Merge the counts back into the original data frame
clinvar.afr <- clinvar.afr %>%
  left_join(gene_counts, by = "GENE") %>%
  rename(gene.counts.AFR = count)
# cc <- cbind.data.frame(clinvar.afr$ID, clinvar.afr$AF, clinvar.afr$AF_afr, clinvar.afr$AF_afr, clinvar.afr$CLNSIG, clinvar.afr$gene.counts.afr, clinvar.afr$ANN....GENE)
dim(clinvar.afr)
# 223

clinvar.afr <- clinvar.afr[clinvar.afr$gene.counts.AFR >= 2,]
dim(clinvar.afr)
# 14

table(clinvar.afr$from)

write.table(as.data.frame(clinvar.afr$ID), "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/vars.clinvar.missense.loftee.afr.maf0.05.mac3.txt", col.names = F, row.names = F, quote = F)

## Write genesets
AFR <- cbind.data.frame(CHROM=clinvar.afr$CHROM, GENE=clinvar.afr$GENE, SNPID=clinvar.afr$ID)
write.table(AFR, "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/extract-variants-final-AFR.GENE-SNP", sep = "\t", row.names = F, quote = F, col.names = F)



