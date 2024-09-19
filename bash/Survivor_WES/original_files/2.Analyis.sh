## Extracting variants on completely QCed (including GQ, DP and VQSR) data
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife
## This is the QC'ed data for Survivors
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated* .

# Get rare_variants_to_extract.txt from extract_PLP.R
plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated --extract rare_variants_to_extract.txt --max-maf 0.01 --keep-allele-order --make-bed -out all_rare_variants_maf0.01_all_sjlife


plink --bfile all_BCC_rare_variants --recodeA --out all_BCC_rare_variants_recodeA
