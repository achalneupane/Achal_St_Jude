## LD check (based on Kateryna's email on 03/17/2023)
# 2 (window size kb, needs to be more than 2) 1 (step size, it could be 1 or 2 in such a small dataset) and 0.8 (it is r2 threshold), generated  plink.prune.in will contain variants that pass the threshold.

# /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ # (old)
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/haplotype_analysis_v2 # (on 07/05/2023)

plink --bfile ../sjlife_ccss_org_ccss_exp_samples --extract ld_list.txt --make-bed --out ttn_significant_ld


plink --bfile ttn_significant_ld --indep-pairwise 2 1 0.2


# plink --bfile ttn_significant_ld --extract plink.prune.in --keep-allele-order --recodeA --out haplotype_input
plink --bfile ttn_significant_ld --extract plink.prune.in --keep-allele-order --recodeA --out haplotype_input_0.2

## Run convert_raw_to_phasing_v2.R script; then run PHASE; 
# PHASE haplotype_input_edited.txt haplotype_phase.out
# PHASE haplotype_input_edited_0.2.txt haplotype_phase_0.2.out
PHASE haplotype_input_edited_0.2.txt haplotype_phase_0.2.out
# an r^2 threshold of 0.8 is often used in LD pruning to retain only variants that are relatively independent. This means that if the r^2 value between two variants is below 0.8, they are considered to be in low LD and are retained, while pairs of variants with an r^2 value equal to or above 0.8 are considered to be in high LD and are pruned.

# then run extract_haplotypes.py