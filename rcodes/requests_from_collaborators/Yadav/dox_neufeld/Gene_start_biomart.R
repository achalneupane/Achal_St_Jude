library(biomaRt)

# Select Ensembl dataset
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", version = "GRCh38")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define gene list
genes <- c("ADD3", "DUSP5", "RBM20", "PDCD4", "BBIP1", "SHOC2", 
           "ADRA2A", "ACSL5", "VTI1A", "TCF7L2", "NRAP", "CASP7", 
           "NHLRC2", "ADRB1", "AFAP1L2", "ABLIM1", "GFRA1", "PNLIPRP2", 
           "HSPA12A", "SLC18A2", "PDZD8", "RAB11FIP2", "MYCBP2", 
           "SLAIN1", "EDNRB", "GPC6", "TGDS", "GPR180", "ABCC4", 
           "UGGT2", "MBNL2")

genes <- as.character(ACMG$GENE)

# Query Ensembl for start and end positions
results <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), 
                 filters = 'hgnc_symbol', 
                 values = genes, 
                 mart = ensembl)

# View the results
print(results)
results <- results[!grepl("PATCH|HSC", results$chromosome_name),]
results$chromosome_name <- paste0("chr",results$chromosome_name)
View(results[-1])

genes[!genes %in% results$hgnc_symbol]

# hgnc_symbol	chromosome_name	start_position	end_position
# ABCC4	13	95019835	95301475
# ABLIM1	10	114431112	114768061
# ACSL5	10	112374116	112428379
# ADD3	10	109996368	110135565
# ADRA2A	10	111077029	111080907
# ADRB1	10	114043866	114046904
# AFAP1L2	10	114294824	114404756
# BBIP1	10	110898730	110919201
# CASP7	10	113679162	113730907
# DUSP5	10	110497907	110511533
# EDNRB	13	77895481	77975529
# GFRA1	10	116056925	116276803
# GPC6	13	93226807	94408020
# GPR180	13	94601857	94634661
# HSPA12A	10	116671192	116850251
# MBNL2	13	97221434	97394120
# MYCBP2	13	77042474	77327094
# NHLRC2	10	113854661	113917194
# NRAP	10	113588714	113664070
# PDCD4	10	110871795	110900006
# PDZD8	10	117277274	117375440
# PNLIPRP2	10	116620953	116645143
# RAB11FIP2	10	118004916	118046941
# RBM20	10	110644336	110839468
# SHOC2	10	110919367	111017307
# SLAIN1	13	77697687	77764242
# SLC18A2	10	117241093	117279430
# TCF7L2	10	112950247	113167678
# TGDS	13	94574054	94596242
# UGGT2	13	95801580	96053482
# VTI1A	10	112446998	112818744

