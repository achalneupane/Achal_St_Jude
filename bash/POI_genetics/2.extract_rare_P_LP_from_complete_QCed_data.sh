cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/POI_genetics
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated_unique.* .


## Keeping QCed SJLIFE, Control and CCSS expansion
awk '{print $1, $2}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam > sjlife_ccss_samples.txt
awk '{print $1, $2}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/CCSS_exp/chr.ALL.CCSS_exp_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam >> sjlife_ccss_samples.txt
awk '{print $1, $2}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam >> sjlife_ccss_samples.txt

## 1. Extract rare variants from Ke et al
module load plink/1.90b
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated_unique --extract Ke_et_al_POI_rare_variants_to_extract.txt --keep sjlife_ccss_samples.txt --make-bed --keep-allele-order --out Ke_et_al_POI_rare_variants_ALL
plink --bfile Ke_et_al_POI_rare_variants_ALL --recodeA --out Ke_et_al_POI_rare_variants_ALL_recodeA


## 2. Extract rare variants from Shekari et al
module load plink/1.90b
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated_unique --extract Shekari_et_al_POI_rare_variants_to_extract.txt --keep sjlife_ccss_samples.txt --make-bed --keep-allele-order --out Shekari_et_al_POI_rare_variants_ALL
plink --bfile Shekari_et_al_POI_rare_variants_ALL --recodeA --out Shekari_et_al_POI_rare_variants_ALL_recodeA


