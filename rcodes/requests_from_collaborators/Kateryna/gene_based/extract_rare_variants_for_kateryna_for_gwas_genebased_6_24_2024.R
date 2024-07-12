# Main folder: /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
clinvar$ID <- sub(";.*", "", clinvar$ID)
clinvar <- clinvar[!duplicated(clinvar$ID),]
# cc <- cbind.data.frame(clinvar$ID, clinvar$AF, clinvar$AF_nfe, clinvar$AF_afr, clinvar$CLNSIG)

#############
## EUR PLP ##
#############
clinvar.nfe <- clinvar[clinvar$AF < 0.01 & clinvar$AF_nfe < 0.01,]
dim(clinvar.nfe)
# 807
## Read MAC data with 0.01 maf
EUR_MAC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/tmp_EUR_merged_data_plink_MAC_EUR.acount", header = T, stringsAsFactors = F, skip = 'CHROM')
## Selecting variants that have MAC of at least 3
EUR_MAC <- EUR_MAC[EUR_MAC$ALT_CTS >= 3,]
clinvar.nfe <- clinvar.nfe[clinvar.nfe$ID %in% EUR_MAC$ID,] # These are also MAF < 0.01 in our EUR data
dim(clinvar.nfe)
# 67

# Keeping only those with 2 or more variants in each gene 
# Count the occurrences of each gene
gene_counts <- clinvar.nfe %>%
  group_by(`ANN....GENE`) %>%
  summarise(count = n())

# Merge the counts back into the original data frame
clinvar.nfe <- clinvar.nfe %>%
  left_join(gene_counts, by = "ANN....GENE") %>%
  rename(gene.counts.NFE = count)
# cc <- cbind.data.frame(clinvar.nfe$ID, clinvar.nfe$AF, clinvar.nfe$AF_nfe, clinvar.nfe$AF_afr, clinvar.nfe$CLNSIG, clinvar.nfe$gene.counts.NFE, clinvar.nfe$ANN....GENE)
dim(clinvar.nfe)
# 67

clinvar.nfe <- clinvar.nfe[clinvar.nfe$gene.counts.NFE >= 2,]
dim(clinvar.nfe)
# 10

write.table(as.data.frame(clinvar.nfe$ID), "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/vars.clinvar.nfe.maf0.01.mac3.txt", col.names = F, row.names = F, quote = F)

## Write genesets
EUR <- cbind.data.frame(CHROM=clinvar.nfe$CHROM, GENE=clinvar.nfe$ANN....GENE, SNPID=clinvar.nfe$ID)
write.table(EUR, "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/extract-variants-final-EUR.GENE-SNP", sep = "\t", row.names = F, quote = F, col.names = F)



#############
## AFR PLP ##
#############
clinvar.afr <- clinvar[clinvar$AF < 0.01 & clinvar$AF_nfe < 0.01,]
dim(clinvar.afr)
# 807
## Read MAC data with 0.01 maf
AFR_MAC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/tmp_AFR_merged_data_plink_MAC_AFR.acount", header = T, stringsAsFactors = F, skip = 'CHROM')
## Selecting variants that have MAC of at least 3
AFR_MAC <- AFR_MAC[AFR_MAC$ALT_CTS >= 3,]
clinvar.afr <- clinvar.afr[clinvar.afr$ID %in% AFR_MAC$ID,] # These are also MAF < 0.01 in our AFR data
dim(clinvar.afr)
# 3

# Keeping only those with 2 or more variants in each gene 
# Count the occurrences of each gene
gene_counts <- clinvar.afr %>%
  group_by(`ANN....GENE`) %>%
  summarise(count = n())

# Merge the counts back into the original data frame
clinvar.afr <- clinvar.afr %>%
  left_join(gene_counts, by = "ANN....GENE") %>%
  rename(gene.counts.NFE = count)
# cc <- cbind.data.frame(clinvar.afr$ID, clinvar.afr$AF, clinvar.afr$AF_nfe, clinvar.afr$AF_afr, clinvar.afr$CLNSIG, clinvar.afr$gene.counts.NFE, clinvar.afr$ANN....GENE)
dim(clinvar.afr)
# 3

clinvar.afr <- clinvar.afr[clinvar.afr$gene.counts.NFE >= 2,]
dim(clinvar.afr)
# 0

write.table(as.data.frame(clinvar.afr$ID), "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/gwas/rare_variants/geneBased_June_24_2024/vars.clinvar.afr.maf0.01.mac3.txt", col.names = F, row.names = F, quote = F)




