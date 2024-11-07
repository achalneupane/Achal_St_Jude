cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ACMG
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated_unique.* .

## Keeping QCed SJLIFE, Control and CCSS expansion
awk '{print $1, $2}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam > sjlife_ccss_samples.txt
awk '{print $1, $2}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/CCSS_exp/chr.ALL.CCSS_exp_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam >> sjlife_ccss_samples.txt
awk '{print $1, $2}' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam >> sjlife_ccss_samples.txt

## Extract rare variants
module load plink/1.90b
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated_unique --extract ACMG_rare_variants_to_extract.txt --keep sjlife_ccss_samples.txt --make-bed --keep-allele-order --out ACMG_rare_variants_ALL
plink --bfile ACMG_rare_variants_ALL --recodeA --out ACMG_rare_variants_ALL_recodeA



