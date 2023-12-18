Gene_set1 <- c("AIP", "AKT1", "ALK", "ANKRD26", "APC", "ATM", "AXIN2", "BAP1", "BARD1", "BLM", "BMPR1A", "BRAF", "BRCA1", "BRCA2", "BRIP1", "BUB1B", "CBL", "CDC73", "CDH1", "CDK4", "CDKN1B", "CDKN1C", "CDKN2A", "CEBPA", "CHEK2", "CREBBP", "CTC1", "CYLD", "DDB2", "DICER1", "DIS3L2", "DKC1", "EGFR", "EP300", "EPCAM", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC6", "ERCC8", "ETV6", "EZH2", "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL", "FANCM", "FH", "FLCN", "GATA2", "GPC3", "HRAS", "KIT", "KRAS", "LIG4", "LZTR1", "MAP2K1", "MAP2K2", "MAX", "MEN1", "MET", "MITF", "MLH1", "MRE11A", "MSH2", "MSH6", "MUTYH", "NBN", "NF1", "NF2", "NHP2", "NRAS", "NSD1", "PALB2", "PDGFRA", "PDGFRB", "PHOX2B", "PIK3CA", "PMS2", "POLD1", "POLE", "POLH", "POT1", "PRKAR1A", "PTCH1", "PTEN", "PTPN11", "RAD51C", "RAD51D", "RAF1", "RB1", "RECQL4", "RET", "RIT1", "RPL11", "RPL35A", "RPL5", "RPS10", "RPS17", "RPS19", "RPS24", "RPS26", "RPS7", "RTEL1", "RUNX1", "SBDS", "SCG5/GREM1", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SH2D1A", "SHOC2", "SLX4", "SMAD4", "SMARCA4", "SMARCB1", "SMARCE1", "SOS1", "SRP72", "STK11", "SUFU", "TERC", "TERT", "TINF2", "TMEM127", "TP53", "TSC1", "TSC2", "VHL", "WAS", "WRAP53", "WRN", "WT1", "XPA", "XPC")
# 73 ACMG genes
ACMG <- c("APC", "RET", "BRCA1", "BRCA2", "PALB2", "SDHD", "SDHAF2", "SDHC", "SDHB", "MAX", "TMEM127", "BMPR1A", "SMAD4", "TP53", "MLH1", "MSH2", "MSH6", "PMS2", "MEN1", "MUTYH", "NF2", "STK11", "PTEN", "RB1", "TSC1", "TSC2", "VHL", "WT1", "FBN1", "TGFBR1", "TGFBR2", "SMAD3", "ACTA2", "MYH11", "PKP2", "DSP", "DSC2", "TMEM43", "DSG2", "RYR2", "CASQ2", "TRDN", "TNNT2", "LMNA", "FLNC", "TTN", "COL3A1", "LDLR", "APOB", "PCSK9", "MYH7", "MYBPC3", "TNNI3", "TPM1", "MYL3", "ACTC1", "PRKAG2", "MYL2", "KCNQ1", "KCNH2", "SCN5A", "BTD", "GLA", "OTC", "GAA", "HFE", "ACVRL1", "ENG", "RYR1", "CACNA1S", "HNF1A", "RPE65", "ATP7B")

unique.genes <- unique(Gene_set1, ACMG)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(Gene_set1 = Gene_set1, ACMG = ACMG),
  category.names = c("Gene_set1", "ACMG"),
  filename = NULL,
  output = TRUE
)

x <- list(Gene_set1=Gene_set1 , ACMG=ACMG)
ggvenn(x, show_elements = F, label_sep = "\n", fill_color = brewer.pal(name="Set2",n=3), text_size= 8, show_percentage = F)
ggvenn(x, show_elements = T, label_sep = "\n", fill_color = brewer.pal(name="Set2",n=3), text_size= 2, auto_scale = TRUE)


Gene_set1[Gene_set1%in% ACMG]
# [1] "APC"     "BMPR1A"  "BRCA1"   "BRCA2"   "MAX"     "MEN1"    "MLH1"    "MSH2"    "MSH6"    "MUTYH"   "NF2"     "PALB2"  
# [13] "PMS2"    "PTEN"    "RB1"     "RET"     "SDHAF2"  "SDHB"    "SDHC"    "SDHD"    "SMAD4"   "STK11"   "TMEM127" "TP53"   
# [25] "TSC1"    "TSC2"    "VHL"     "WT1"  