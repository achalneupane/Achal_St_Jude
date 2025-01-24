## Kateryn's GWAS project
## Top SNP: chr16:25611595:C:T
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WGS_Northwestern_Joint_call_153_samples/Kateryna_gwas
# module load plink/1.90b
# plink --chr 16 --from-bp 25611595 --make-bed --out chr16_25611595_C_T --to-bp 25611595 --double-id --vcf-half-call m --keep-allele-order --vcf ../CAB_5428.haplotype.recalibrated.decomposed_chr16.vcf.gz
# plink --bfile chr16_25611595_C_T --recodeA --keep-allele-order --out chr16_25611595_C_T_recodeA 

RNAseq <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/eQTL_RNAseq_NW_12_6_2024/SAPKO-842478-UNSTRANDED_RSEM_gene_TPM.2025-01-07_19-54-28.txt", header = T)
raw <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/WGS_Northwestern_Joint_call_153_samples/Kateryna_gwas/chr16_25611595_C_T_recodeA.raw")

rnaseq.samples <- as.data.frame(colnames(RNAseq))
rnaseq.samples$IID2 <- sub("_.*", "", rnaseq.samples$`colnames(RNAseq)`)

colnames(RNAseq)
# [1] "geneID"          "geneSymbol"      "bioType"         "annotationLevel" "CE01_0"          "CE01_1"          "CE01_3"          "JR06_0"          "JR06_1"         
# [10] "JR06_3"          "JR12c2_0"        "JR12c2_1"        "JR12c2_3"        "JR20c1_0"        "JR20c1_1"        "JR20c1_3"        "JR53c8_0"        "JR53c8_1"       
# [19] "JR53c8_3"        "JR54c1_0"        "JR54c1_1"        "JR54c1_3"        "JR63_0"          "JR63_1"          "JR63_3"          "JR71_0"          "JR71_1"         
# [28] "JR71_3"          "KG13_13_0"       "KG13_13_1"       "KG13_13_3"       "KG20_7_0"        "KG20_7_1"        "KG20_7_3"        "KG23c12_0"       "KG23c12_1"      
# [37] "KG23c12_3"       "KG33_0"          "KG33_1"          "KG33_3"          "KG34C1_0"        "KG34C1_1"        "KG34C1_3"        "KG37_8_0"        "KG37_8_1"       
# [46] "KG37_8_3"        "PT01_0"          "PT01_1"          "PT01_3"          "SC01_c15_0"      "SC01_c15_1"      "SC01_c15_3"      "SC02_6_0"        "SC02_6_1"       
# [55] "SC02_6_3"        "SC04c2_0"        "SC04c2_1"        "SC04c2_3"        "SC06_2_0"        "SC06_2_1"        "SC06_2_3"        "SC12_c5_0"       "SC12_c5_1"      
# [64] "SC12_c5_3"       "SC13_16_0"       "SC13_16_1"       "SC13_16_3"       "SC14_5_0"        "SC14_5_1"        "SC14_5_3"        "SC16_4_0"        "SC16_4_1"       
# [73] "SC16_4_3"        "SC18_2_0"        "SC18_2_1"        "SC18_2_3"        "SC23_0"          "SC23_1"          "SC23_3"          "SC281A_0"        "SC281A_1"       
# [82] "SC281A_3"        "SC29_11_0"       "SC29_11_1"       "SC29_11_3"       "SC34_0"          "SC34_1"          "SC34_3"          "SC35_0"          "SC35_1"         
# [91] "SC35_3"          "SC41c11_0"       "SC41c11_1"       "SC41c11_3"      

# SC281A_0 is SC28

unique(raw$IID)
unique(raw$IID)
# [1] "442332"       "442333c13_5"  "442333c8"     "CE01_C5_P26"  "CE11-5"       "CE12-1"       "CE16-1"       "CE17"         "CE20-2"       "CE21-14"      "CE22-1"      
# [12] "CE4-7"        "CE5"          "CE6"          "CE7"          "CE8"          "CE9"          "GW10-9"       "GW115-8"      "GW124-9"      "GW129-24"     "GW132-2B"    
# [23] "GW159-11"     "GW165-11"     "GW167-6"      "GW168-1"      "GW169-31"     "GW2-6"        "GW21-6"       "GW22-36"      "GW28-12"      "GW29-3"       "GW3-7"       
# [34] "GW30-14"      "GW4-12"       "GW53-5"       "GW64-1"       "HP20"         "HP21"         "JR01-C33"     "JR06-C5-p21"  "JR07-C11-p30" "JR08-C2-p20"  "JR10-C1-p26" 
# [45] "JR11-C8-p24"  "JR12-C2-p27"  "JR17-C1-p40"  "JR18_C7_P20"  "JR19-C3-p22"  "JR20_C1_P24"  "JR23_C6_P18"  "JR25-C5-p24"  "JR28-C5-p21"  "JR33-C9-p19"  "JR37-C2--p21"
# [56] "JR44-C4-p21"  "JR45_C15_P20" "JR51-C2-p22"  "JR52-C6-p25"  "JR53_C8_P24"  "JR54_C5_P23"  "JR55_C11_P9"  "JR56_C5_P12"  "JR57_C2_P33"  "JR58-C7-p44"  "JR58-C7-p48" 
# [67] "JR59_C1_P20"  "JR60-C11-p20" "JR63_C21_P20" "JR64_C5_P26"  "JR65_C3_P16"  "JR66_C3_P24"  "JR68_C8_P42"  "JR69-C3-p19"  "JR71-C17-p25" "JW1-12"       "JW10-1"      
# [78] "JW11-6"       "JW12-17"      "JW13-2"       "JW2-2"        "JW3-4"        "JW4-5"        "JW5-5"        "JW6"          "JW7-7"        "JW8-6"        "JW9-12"      
# [89] "KG01-C5-p26"  "KG03-C13-p22" "KG07-C6-p22"  "KG09-C8-p21"  "KG10-c9-p23"  "KG10_C9_P23"  "KG11-C6-p21"  "KG12-C3-p18"  "KG13_C13_P27" "KG19_C6_P26"  "KG20-C9-p21" 
# [100] "KG23_C12_P19" "KG24_C15_P22" "KG25-C7-p26"  "KG27-C6-p18"  "KG28-C5-p17"  "KG30-C10-p34" "KG32-C6-p23"  "KG33_C2_P26"  "KG34_C1_P20"  "KG36-C10-p20" "KG37-C9-p15" 
# [111] "KG38-C12-p14" "KG39_C6_P20"  "KG41-c7-p13"  "KG42_C15_P24" "PT01-C8-p20"  "SC01-C15-p20" "SC02-C6-P22"  "SC03-C1-p20"  "SC04-C2-p21"  "SC05-C7-p27"  "SC08-C1-p18" 
# [122] "SC09-C5-p11"  "SC10-C14-p12" "SC11-C3-p12"  "SC12-C5-p24"  "SC13-C16-p21" "SC14-C5-p19"  "SC15-C2-p15"  "SC16_C4_P20"  "SC17-C15-p13" "SC18-C10-p16" "SC19-C2-p10" 
# [133] "SC20_C12_P24" "SC21--C3-p13" "SC22-C01-P9"  "SC23-C2-p8"   "SC24-C8-p9"   "SC25-C5-p7"   "SC26-C8-p9"   "SC27-C5-p6"   "SC28_C1_P52"  "SC29_C11_P18" "SC30_C5_P22" 
# [144] "SC31_C1_P15"  "SC32_C12_P36" "SC34_C8_P24"  "SC35"         "SC36_C5_P21"  "SC37_C6_P34"  "SC38_C2_P14"  "SC39_C11_P44" "SC41_C11_P41" "SC6-C2-p28" 

raw$IID2 <- sub("_.*", "", raw$IID)
raw$IID2 <- sub("-.*", "", raw$IID2)

unique(rnaseq.samples$IID2)[!unique(rnaseq.samples$IID2) %in% raw$IID2]

rnaseq.samples$IID2[grepl("JR", rnaseq.samples$IID2)] <- gsub("c1|c2|c8|c12|C1", "", rnaseq.samples$IID2[grepl("JR", rnaseq.samples$IID2)])
rnaseq.samples$IID2[grepl("KG", rnaseq.samples$IID2)] <- gsub("c1|c2|c8|c12|C1", "", rnaseq.samples$IID2[grepl("JR", rnaseq.samples$IID2)])
