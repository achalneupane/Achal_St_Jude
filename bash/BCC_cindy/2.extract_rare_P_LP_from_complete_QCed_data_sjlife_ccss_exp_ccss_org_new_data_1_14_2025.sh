#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
## Extract rare P/LP from QCed data which did go through GQ, DP and VQSR filter
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis2
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/plink_all/chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated.* .
module load plink/1.90b
# ## Extracting variants on GQ, DP and VQSR data only
# plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated --extract rare_variants_to_extract.txt --keep-allele-order --make-bed -out all_BCC_rare_variants_VQSR
# plink --bfile all_BCC_rare_variants_VQSR --recodeA --out all_BCC_rare_variants_VQSR_recodeA


## Extracting variants on completely QCed (including GQ, DP and VQSR) data
plink --bfile chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated --extract rare_variants_to_extract.txt --keep-allele-order --make-bed -out all_BCC_rare_variants
# 8030 people (0 males, 0 females, 8030 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to all_BCC_rare_variants.nosex .
# --extract: 29711 variants remaining.

plink --bfile all_BCC_rare_variants --recodeA --out all_BCC_rare_variants_recodeA



1. For the CSG60 dataset, we found 43/60 genes, with a total of 285 variants in WES data, utilizing ClinVar annotations. Out of the 31 reported genes (Clinvar set), 29 were found. Among these, 67 variants from these genes matched with 148 variants reported by Kim et al in the ClinVar set (21 of which were not considered pathogenic/likely pathogenic (P/LP) in the latest ClinVar update).
2. For the CSG172 dataset, we found 84/172 genes, with a total of 473 variants in WES data. Out of the 46 reported genes (Clinvar set), 41 were found. Among these, 64 variants matched with 188 variants reported by Kim et al in the ClinVar set.


grep 'NA' csg.60.vars.unique.txt | cut -f8 | while read -r line; do
  # echo "Doing ${line}"
  grep "${line}" /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated_unique.bim
done



CSG60 variants present in WES but not P/LP (all uncertain significance)
# 2       chr2:47800130:C:*;rs587782805       0       47800130        *       C
# 17      chr17:31169999:T:G      0       31169999        G       T
# 17      chr17:31227289:G:T      0       31227289        T       G
# 9       chr9:95485875:C:T;rs368869806       0       95485875        T       C
# 9       chr9:37002705:C:T;rs398123063       0       37002705        T       C
# 13      chr13:48360062:T:A      0       48360062        A       T
# 13      chr13:48303931:C:T;rs1952051704      0       48303931        T       C
# 13      chr13:48465238:C:T;rs137853293      0       48465238        T       C
# 13      chr13:48303991:C:CGGAACCCCCGGCA 0       48303991        CGGAACCCCCGGCA  C
# 13      chr13:48456349:G:T      0       48456349        T       G
# 13      chr13:48349000:G:A;rs2138093650      0       48349000        A       G
# 10      chr10:43119576:G:A;rs1318733775      0       43119576        A       G
# 1       chr1:17028712:T:C;rs1310341038       0       17028712        C       T
# 1       chr1:17033077:C:T;rs570278423       0       17033077        T       C
# 1       chr1:241512001:G:C;rs199822819      0       241512001       C       G
# 10      chr10:102509250:G:A;rs2135620546     0       102509250       A       G
# 17      chr17:7675058:C:T;rs150607408       0       7675058 T       C
# 11      chr11:32428554:G:A      0       32428554        A       G


grep 'NA' csg.60.vars.unique.txt | cut -f8 | while read -r line; do
  # echo "Doing ${line}"
  grep "${line}" /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/all.bim
done

CSG60 variants present in WGS (all uncertain significance)
# 5       chr5:112834948:T:G      0       112834948       G       T **
# 1       chr1:241512001:G:C      0       241512001       C       G
# 2       chr2:47800130:C:T       0       47800130        T       C
# 17      chr17:31169999:T:G      0       31169999        G       T
# 17      chr17:31259060:T:C      0       31259060        C       T
# 17      chr17:31227289:G:T      0       31227289        T       G
# 9       chr9:37002705:C:T       0       37002705        T       C
# 9       chr9:95485875:C:T       0       95485875        T       C
# 10      chr10:87863635:C:T      0       87863635        T       C **
# 13      chr13:48360062:T:A      0       48360062        A       T
# 13      chr13:48303931:C:T      0       48303931        T       C
# 13      chr13:48465238:C:T      0       48465238        T       C
# 13      chr13:48456349:G:T      0       48456349        T       G
# 13      chr13:48349000:G:A      0       48349000        A       G
# 13      chr13:48303716:G:T      0       48303716        T       G **
# 10      chr10:43119576:G:A      0       43119576        A       G
# 1       chr1:17028712:T:C       0       17028712        C       T
# 1       chr1:17033077:C:T       0       17033077        T       C
# 1       chr1:170330771:T:C      0       170330771       C       T **
# 11      chr11:112093129:C:T     0       112093129       T       C **
# 11      chr11:112093129:C:<*:DEL>       0       112093129       <*:DEL> C **
# 10      chr10:102509250:G:A     0       102509250       A       G
# 17      chr17:76761531:G:A      0       76761531        A       G **
# 17      chr17:7675058:C:T       0       7675058 T       C
# 17      chr17:76750589:C:T      0       76750589        T       C
# 17      chr17:7675070:C:T       0       7675070 T       C
# 11      chr11:32428554:G:A      0       32428554        A       G



grep 'NA' kim_ST1.csg172.vars.unique.txt | cut -f10 | while read -r line; do
  # echo "Doing ${line}"
  grep -w "${line}" /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim
done


CSG172 variants present in WES but not P/LP (all uncertain significance)
1       chr1:241497927:A:ATTT;rs367543046   0       241497927       ATTT    A ## conflicting
1       chr1:241512001:G:C;rs199822819       0       241512001       C       G ## conflicting
2       chr2:136114925:C:T;rs147214773      0       136114925       T       C ## conflicting
3       chr3:10149921:C:T;rs28940298       0       10149921        T       C ## conflicting
3       chr3:48583584:G:A;rs79378857       0       48583584        A       G ## conflicting
3       chr3:69959395:A:G;rs368915509       0       69959395        G       A
5       chr5:1294282:C:T;rs121918661        0       1294282 T       C ## conflicting
5       chr5:112780865:C:G;rs141576417     0       112780865       G       C ## conflicting
14      chr14:81062179:C:T;rs142063461      0       81062179        T       C ## conflicting
14      chr14:81143695:G:A;rs121908866      0       81143695        A       G ## conflicting
17      chr17:58709883:A:G;rs199886026      0       58709883        G       A ## conflicting
17      chr17:7675070:C:T;rs397514495       0       7675070 T       C ## conflicting
17      chr17:61716051:G:A;59793412       0       61716051        A       G ## conflicting
22      chr22:28695868:AG:A;rs555607708     0       28695868        A       AG ## conflicting
22      chr22:28711986:C:T;rs121908702      0       28711986        T       C ## Uncertain significance
22      chr22:28725242:C:*;      0       28725242        *       C ## conflicting
22      chr22:28725242:C:T;rs121908698      0       28725242        T       C ## conflicting



plink --bfile all_BCC_rare_variants --keep IID_fam_all_EUR_from_subdf.txt --make-bed --out all_BCC_rare_variants_EUR

plink --bfile all_BCC_rare_variants_EUR --freq --out all_BCC_rare_variants_EUR_freq



##################################
## 08/15/2024                   ##
## Discrepancy check for CSG172 ##
##################################

From Kim et al: file:///Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/rare_variants/kim_et_al/pkab007.pdf
Title: Frequency of Pathogenic Germline Variants in Cancer-Susceptibility Genes in the Childhood Cancer Survivor Study

"Methods: Exome sequencing was conducted on germline DNA from 5451 pediatric
cancer survivors (cases who survived 5 years from diagnosis; n ¼ 5105 European) and 597 European cancer-free adults (controls). Analyses focused on comparing the frequency of rare P/LP variants in 237 cancer-susceptibility genes and a subset of
60 autosomal dominant high-to-moderate penetrance genes, for both case-case and case-control comparisons."

"We focused our analyses on 172 autosomal dominant genes, including both autosomal dominant and autosomal recessive genes
and those of unknown inheritance (hereafter, CSG-172), and a subset of 60 autosomal dominant genes with high-to-moderate
penetrance (hereafter, CSG-60) previously reported"

"The analysis was restricted to variants with a minor allele frequency of less than 1% in cases or controls drawn from all ethnic subgroups reported in Exome Aggregation Consortium
(excluding The Cancer Genome Atlas data)"

"Variants were categorized as pathogenic (P), likely pathogenic (LP), variant of uncertain significance, likely benign, or benign using a hierarchical classification system based on ClinVar, the Human
Gene Mutation Database, and InterVar" >>>>> In our case, we are considering Clinvar only

CSG60: Percent P/LP (95% CI) 4.1 (3.6 to 4.7) 1.3 (0.4 to 2.3)
CSG172: Percent P/LP (95% CI) 11.9 (11.0 to 12.8) 7.7 (5.6 to 9.8)   

#################################
## Discrepancy check for CSG60 ##
#################################
Title: "Genetic Risk for Subsequent Neoplasms Among Long-Term Survivors of Childhood Cancer"
"CSG60: Pathogenic/likely pathogenic (P/LP) germline mutations in cancer predisposition genes have been shown to be present in 8.5% of children with newly diagnosed
cancer".
"Participants were 3,006 survivors (53% male; median age, 35.8 years [range, 7.1 to 69.8 years]; 56%
received radiotherapy), 1,120 SNs were diagnosed among 439 survivors (14.6%), and 175 P/LP
mutations were identified in 5.8% (95% CI, 5.0% to 6.7%) of survivors."
Paper: Z:\ResearchHome\Groups\sapkogrp\projects\Genomics\common\BCC\files_shared_by_cindy\rare_variants\CSG60_zhaoming


"By combing the two large cohorts of long-term childhood 
cancer survivors with germline genome or exome 
sequencing data, we comprehensively examined the 
prevalence and spectrum of cancer predisposing variants 
by specific cancer diagnosis and cancer predisposition 
gene. The 6·7% prevalence of cancer predisposing 
variants among the 4402 participants in SJLIFE is slightly 
higher than the 5·8% prevalence that was previously 
reported on the basis of a subset of 3006 participants in SJLIFE, and is substantially higher than the 4·3% 
prevalence in CCSS, primarily because of the different 
composition of childhood cancers in these two cohorts."
Paper: Z:\ResearchHome\Groups\sapkogrp\projects\Genomics\common\BCC\files_shared_by_cindy\rare_variants\Zhaoming_lancet_oncology



