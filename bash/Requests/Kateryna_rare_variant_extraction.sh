cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
## SnpEff annotation (Clinvar and LOF variants)
head -1 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt > LoF_variants.txt
cat $(ls *.dbSNP155-FIELDS-simple.txt| sort -V)| egrep 'frameshift_variant|start_lost|stop_gained|splice|splice_region_variant' >> LoF_variants.txt
awk '{ print $3 }' LoF_variants.txt| cut -d";" -f1 > LoF_variants_ID.txt



## R code
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/")
## read annotated SJLIFE annotated VCF 
# Loop over all chromosomes
chromosomes <- 1:22


## Now extract variants of ClinVar and MetaSVM significance
FINAL.VCF <- {}

capture.output (for( i in 1:length(chromosomes)){
print(paste0("Doing chromosome ", chromosomes[i]))
  
VCF <- fread(paste0("MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr", chromosomes[i], ".preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt"))

#####################################
## Clinvar based P/LP VCF variants ##
#####################################
print(as.data.frame(table(VCF$CLNSIG)))
# Var1     Freq
# 1                                                          . 19388477
# 2                                                    Affects       51
# 3                                                     Benign   117561
# 4                                       Benign/Likely_benign    15102
# 5                                Benign/Likely_benign|_other        7
# 6                                              Benign|_other       14
# 7               Conflicting_interpretations_of_pathogenicity    19244
# 8  Conflicting_interpretations_of_pathogenicity|_association        6
# 9  Conflicting_interpretations_of_pathogenicity|_risk_factor        6
# 10                                             Likely_benign    62840
# 11                                         Likely_pathogenic      593**
# 12                                                Pathogenic     2028 **
# 13                              Pathogenic/Likely_pathogenic      875**
# 14                                   Pathogenic|_risk_factor        9**
# 15                                    Uncertain_significance    47896
# 16                                               association        7
# 17                                              not_provided      535
# 18                                                protective        3
# 19                                   protective|_risk_factor        5
# 20                                               risk_factor       22

## Wanted Clinvar patterns
WANTED.types.clinvar <- c("^Pathogenic/Likely_pathogenic$|^Likely_pathogenic$|^Likely_pathogenic/Pathogenic$|^Pathogenic$|Pathogenic\\|_risk_factor|^Pathogenic\\|")
print(paste0("Total P or LP vars from clinvar: ", sum(grepl(WANTED.types.clinvar, VCF$CLNSIG, ignore.case = T))))

## Wanted MetaSVM patterns D (available patterns: D= Deleterious; T= Tolerated)
wanted.types.MetaSVM <- c("D")
print(paste0("Total Deleterious vars from MetaSVM: ", sum(grepl(wanted.types.MetaSVM, VCF$dbNSFP_MetaSVM_pred, ignore.case = T))))

VCF.clinvar <- VCF[grepl(WANTED.types.clinvar, VCF$CLNSIG, ignore.case = T),]
VCF.clinvar$PRED_TYPE <- "Clinvar"
VCF.MetaSVM <- VCF[grepl(wanted.types.MetaSVM, VCF$dbNSFP_MetaSVM_pred, ignore.case = T),]
VCF.MetaSVM$PRED_TYPE="MetaSVM"
VCF <- rbind.data.frame(VCF.clinvar, VCF.MetaSVM)
FINAL.VCF <- rbind.data.frame(FINAL.VCF, VCF)
}, file = "SNPEFF_clinvar_metaSVM_from_R_filtering_process.log")


# # First remove any leading and trailing spaces
FINAL.VCF$CHROM <- trimws(FINAL.VCF$CHROM, which = "both")
FINAL.VCF$POS <- trimws(FINAL.VCF$POS, which = "both")

FINAL.VCF$KEY <- paste(FINAL.VCF$CHROM, FINAL.VCF$POS, FINAL.VCF$REF, FINAL.VCF$ALT, sep =":")

# FINAL.VCF.unique <- FINAL.VCF[duplicated(FINAL.VCF$KEY),]

#############
## Clinvar ##
#############
CLINVAR <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "Clinvar",] ## You need this
table(CLINVAR$CLNSIG)
nrow(CLINVAR)
# 60892

#############
## MetaSVM ##
#############
MetaSVM <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "MetaSVM",]

CLINVAR.unique <- CLINVAR[-which(duplicated(CLINVAR$KEY)),]
table(CLINVAR.unique$CLNSIG)

MetaSVM.unique <- MetaSVM[!duplicated(MetaSVM$KEY),]

nrow(CLINVAR.unique)
# 4705
nrow(MetaSVM.unique)
# 83479

######################
## Loss of Function ##
######################
LoF <- read.delim("LoF_variants.txt",  sep = "\t", header = T, check.names = F)
LoF$PRED_TYPE <- "LoF"
LoF$KEY <- paste(LoF$CHROM, LoF$POS, LoF$REF, LoF$ALT, sep = ":")
LoF.unique <- LoF[!duplicated(LoF$KEY),]
nrow(LoF.unique)
# 316386

splice.region <- LoF.unique[grepl("splice_region_variant", LoF.unique$`ANN[*].EFFECT`),]
dim(splice.region)
# 242541     25
non.splice.region <- LoF.unique[!grepl("splice_region_variant", LoF.unique$`ANN[*].EFFECT`),]
dim(non.splice.region)
# 73845    25

NO.splice.region <- splice.region[grepl("frameshift|stop|donor|acceptor", splice.region$`ANN[*].EFFECT`),]
dim(NO.splice.region)
# 4421   25

LOF.keep <- rbind.data.frame(non.splice.region, NO.splice.region)
dim(LOF.keep)

Predicted.vars.inVCF.Unique <-  rbind.data.frame(CLINVAR.unique, LOF.keep)





# Then, get the gnomAD MAF for EUR and ALL samples.
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar
# For chr17: ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr17.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt

## Extract gnomAD with maf < 0.01 lines
head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt | sed 's/\t/\n/g' | nl
# head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01.txt

for i in {1..22}; do
CHROM="chr${i}"
echo "Doing ${CHROM}"
head -1 /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt| awk '{print $178"\t"$179"\t"$181"\t"$182"\t"$38"\t"$39"\t"$44}' > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01_${CHROM}.txt	
awk -F'\t' '{if($38 != "." && $38 < 0.01) print $178"\t"$179"\t"$181"\t"$182"\t"$38"\t"$39"\t"$44}' "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/ANNOVAR_MERGED.SJLIFE.1.2.GATKv3.4.VQSR.${CHROM}.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt" >> /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01_${CHROM}.txt
done




gnomAD_genome_ALL gnomAD_exome_NFE gnomAD_genome_AFR 

df <- df[df$gnomAD_genome_ALL < 1% & df$gnomAD_exome_NFE < 1%,]




# Extract haplotype
cd  /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/haplotype
cut -d' ' -f1,2 pheno/sjlife_afr_dox_only_pcs.pheno > haplotype/AFRsamples
plink --vcf ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz --make-bed --keep EURsamples --extract significant --out haplotype_input_eur
plink --vcf ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz --make-bed --keep AFRsamples --extract significant --out haplotype_input_afr











# Haplotype
plink --bfile haplotype_input_afr --keep-allele-order --recodeA --out haplotype_input_afr_phase_raw
plink --bfile haplotype_input_eur --keep-allele-order --recodeA --out haplotype_input_afr_phase_eur


PHASE haplotype_input_edited_0.2.txt haplotype_phase_0.2.out



## Extract clinvar
cut -d$'\t' -f22 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD.gnomAD_clinvar_12_10_2023-clinvar_12_10_2023_FIELDS-simple.txt| sort -V | uniq
awk -F'\t' 'NR==1 || ($29 ~ /Pathogenic|Likely_pathogenic/) && $29 !~ /Conflicting|Benign/' MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD.gnomAD_clinvar_12_10_2023-clinvar_12_10_2023_FIELDS-simple.txt > all_plp_kateryna_dyslipidemia_new_clinvar.txt

# extract MAF < 1%; 54th column is AF
awk 'NR == 1 || $54 < 0.01' all_plp_kateryna_dyslipidemia_new_clinvar.txt > all_plp_kateryna_dyslipidemia_new_clinvar_0.01_AF.txt

## Now extract genes
<file_start_end.txt>
chr9	10105	14105	ABC
chr9	18105	24105	DEF

while read -r CHR START END GENE; do
  # Extract the header from the input file
  awk 'NR==1 {print}' all_plp_kateryna_dyslipidemia_new_clinvar_0.01_AF.txt > "result_${GENE}_new_clinvar.txt"
  # Filter rows based on conditions and append to the result file
  awk -v chr="$CHR" -v start="$START" -v end="$END" '$1 == chr && $2 >= start && $2 <= end' all_plp_kateryna_dyslipidemia_new_clinvar_0.01_AF.txt  >> "result_${GENE}_new_clinvar.txt"
done < file_start_end.txt


cut -d$'\t' -f22 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD.gnomAD_clinvar_12_10_2023-clinvar_12_10_2023_FIELDS-simple.txt| sort -V | uniq -c
# 585615527 .
#     237 Affects
#       5 Affects_association
#       1 Affects_risk_factor
#  198697 Benign
#    2265 Benign,Benign
#      34 Benign,Benign,Benign
#   90505 Benign/Likely_benign
#    1208 Benign/Likely_benign,Benign/Likely_benign
#       7 Benign/Likely_benign_association
#       9 Benign/Likely_benign_drug_response
#       3 Benign/Likely_benign_drug_response_other
#      65 Benign/Likely_benign_other
#      13 Benign/Likely_benign_other_risk_factor
#      43 Benign/Likely_benign_risk_factor
#       3 Benign_confers_sensitivity
#      14 Benign_drug_response
#     161 Benign_other
#       8 Benign_protective
#      11 Benign_risk_factor
#  250614 Conflicting_interpretations_of_pathogenicity
#    2571 Conflicting_interpretations_of_pathogenicity,Conflicting_interpretations_of_pathogenicity
#      41 Conflicting_interpretations_of_pathogenicity_Affects
#      68 Conflicting_interpretations_of_pathogenicity_association
#      13 Conflicting_interpretations_of_pathogenicity_association_risk_factor
#      17 Conflicting_interpretations_of_pathogenicity_drug_response
#     335 Conflicting_interpretations_of_pathogenicity_other
#      24 Conflicting_interpretations_of_pathogenicity_other_risk_factor
#      99 Conflicting_interpretations_of_pathogenicity_risk_factor
#  153356 Likely_benign
#    2775 Likely_benign,Likely_benign
#      18 Likely_benign,Likely_benign,Likely_benign
#       3 Likely_benign_association
#      54 Likely_benign_drug_response_other
#     107 Likely_benign_other
#    8627 Likely_pathogenic
#     138 Likely_pathogenic,Likely_pathogenic
#      10 Likely_pathogenic_risk_factor
#      32 Likely_risk_allele
#   20805 Pathogenic
#     275 Pathogenic,Pathogenic
#   14533 Pathogenic/Likely_pathogenic
#     203 Pathogenic/Likely_pathogenic,Pathogenic/Likely_pathogenic
#      15 Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance
#      30 Pathogenic/Likely_pathogenic_other
#      33 Pathogenic/Likely_pathogenic_risk_factor
#       2 Pathogenic/Pathogenic,_low_penetrance_risk_factor
#      13 Pathogenic_Affects
#      17 Pathogenic_association
#      29 Pathogenic_drug_response
#      18 Pathogenic_other
#      16 Pathogenic_protective
#      40 Pathogenic_risk_factor
#      26 Uncertain_risk_allele
#      10 Uncertain_risk_allele_protective
#       2 Uncertain_risk_allele_risk_factor
# 1043649 Uncertain_significance
#   14156 Uncertain_significance,Uncertain_significance
#     176 Uncertain_significance,Uncertain_significance,Uncertain_significance
#      82 Uncertain_significance/Uncertain_risk_allele
#       8 Uncertain_significance_association
#      18 Uncertain_significance_drug_response
#      12 Uncertain_significance_other
#       8 Uncertain_significance_risk_factor
#     240 association
#       1 association_drug_response
#      12 confers_sensitivity
#      22 dbNSFP_clinvar_clnsig
#     716 drug_response
#       3 drug_response_other
#      22 drug_response_risk_factor
#    2539 not_provided
#      19 not_provided,not_provided
#      72 other
#       6 other,other
#      55 protective
#      12 protective_risk_factor
#     624 risk_factor
