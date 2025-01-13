cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/POI_genetics
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/plink_all/chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated.* .



## 1. Extract rare variants from Ke Shekari et al
module load plink/1.90b
plink --bfile chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated --extract Ke_Shekari_et_al_POI_rare_variants_to_extract.txt --make-bed --keep-allele-order --out Ke_Shekari_et_al_POI_rare_variants_ALL
plink --bfile Ke_Shekari_et_al_POI_rare_variants_ALL --recodeA --out Ke_Shekari_et_al_POI_rare_variants_ALL_recodeA




