## LD check (based on Kateryna's email on 03/17/2023)
# 2 (window size kb, needs to be more than 2) 1 (step size, it could be 1 or 2 in such a small dataset) and 0.8 (it is r2 threshold), generated  plink.prune.in will contain variants that pass the threshold.

# /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ # (old)
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/haplotype_analysis_v2 # (on 07/05/2023)

plink --bfile ../sjlife_ccss_org_ccss_exp_samples_v2 --extract ld_list.txt --make-bed --out ttn_significant_ld


plink --bfile ttn_significant_ld --indep-pairwise 2 1 0.2


# plink --bfile ttn_significant_ld --extract plink.prune.in --keep-allele-order --recodeA --out haplotype_input
plink --bfile ttn_significant_ld --extract plink.prune.in --keep-allele-order --recodeA --out haplotype_input_0.2

## Run convert_raw_to_phasing_v2.R script; then run PHASE; 
# PHASE haplotype_input_edited.txt haplotype_phase.out
# PHASE haplotype_input_edited_0.2.txt haplotype_phase_0.2.out
PHASE haplotype_input_edited_eur.txt haplotype_input_edited_eur.txt.out

# then run extract_haplotypes.py; Then run association_analysis_with_haplotypes_v2.R

# In haplotype association analysis, the analysis can be performed for variants that are either in low linkage disequilibrium (LD) or high linkage disequilibrium, depending on the specific research question and study design. The choice of including variants in low or high LD depends on the goals of the analysis and the specific hypotheses being tested.

# Here are two scenarios to consider:

# Including Variants in Low LD:

# If the goal is to identify haplotypes that are directly associated with a phenotype or disease outcome, it might be more appropriate to include variants in low LD.
# Including variants in low LD ensures that each variant is providing independent or unique genetic information, reducing the possibility of redundancy in the analysis.
# This approach allows for the identification of specific haplotypes that might have a direct impact on the phenotype of interest.
# Including Variants in High LD:

# If the focus is on capturing the overall LD structure or haplotype blocks in the genomic region, including variants in high LD can be beneficial.
# Haplotypes in high LD tend to co-occur and are inherited together as a unit, representing a broader pattern of genetic variation.
# Including variants in high LD can provide a more comprehensive picture of the LD structure, allowing for the identification of haplotype blocks or regions that are associated with the phenotype or disease outcome.
# In practice, the choice of including variants in low or high LD depends on the research question, the specific genetic region under investigation, the availability of relevant prior knowledge, and the study design. Researchers need to carefully consider the biological context and hypotheses being tested when deciding which LD patterns to include in the haplotype association analysis.