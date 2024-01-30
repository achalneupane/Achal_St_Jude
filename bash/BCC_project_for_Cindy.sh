# Dear Achal,

# I am writing to ask you for your help. Do you think you might be able to help us put the following together over the next 2 weeks to meet the mid-Feb ISLCCC abstract deadline?

# Basically, we would like to get P/LP variant carrier variables in SJLIFE (WES/WGS), CCSS Expansion (WES/WGS), and CCSS Original (WES) considering the following 5 gene lists below, using any of the 3 rare variant masks that were described in our POI concept proposal (SnpEff; LOFTEE; ClinVar - I think if we have to pick one, we should pick the strictest definition, which is probably ClinVar). Since this is a quick preliminary pass, maybe we could focus on getting at least one cancer susceptibility P/LP variant carrier status variable ([1] or [2] below) and at least one BCC-related variable ([3]-[5] below). Do you think this would be feasible? Please let me know if you have questions or if you would like to meet to discuss.

# Thank you so much,
# Cindy

# Gene lists 

# (1) 60 cancer susceptibility genes - should be Zhaoming's SJCPG60 list. See Kim_ST1.txt, use "CSG_60" (genes marked with "x"). From Kim et al evaluating cancer susceptibility gene P/LP variants in CCSS, link: https://doi.org/10.1093/jncics/pkab007.

# (2) Expanded list of 172 cancer susceptibility genes. Kim_ST1.txt, use "CSG_172" (genes marked with "x"; 172 genes evaluated in CCSS)

# (3) Literature-based list of BCC genes. See "NCI table 3" for a list of ~20 genes (Table 3 in the cancer.gov). These genes are based on what has been reported in the literature re: BCC-associated syndromes.
# Based on Kilgour et al, Choquet et al, and NCI cancer.gov (all 3 are very consistent):
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8345475/#B6-cancers-13-03870
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7259534/
# https://www.cancer.gov/types/skin/hp/skin-genetics-pdq

# (4) BCC gene panel genes. See "BCC_blueprint_panel.txt" for a list of genes tested by Blueprint Genetics panel for BCC:
# https://blueprintgenetics.com/tests/panels/dermatology/hereditary-melanoma-and-skin-cancer-panel/

# (5) ClinVar BCC-related genes. See "clinvar_result_BCC.txt". These are all ClinVar listings for "basal AND cell AND carcinoma" annotated as P or LP and with multiple submitters/no conflicts. This list hasn't been reviewed yet in terms of conditions, but I think we can use the variant/gene list. For the final analysis, I can ask the clinicians to help review the conditions to make sure we don't have anything unrelated.


# Further details in Email from 1/25/2024


## Extract variants for PRS
# load modules
module load bcftools
module load plink/1.90b


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/prs/sjlife
## SJLIFE
THREADS=4
for CHR in {1..22}; do
grep -w chr${CHR} ../all_bed_BCC.bed | sort -V > PRS_vars_chr${CHR}.bed
sed -i "s/\r//g" PRS_vars_chr${CHR}.bed
bcftools view -Oz /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2//MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.PASS.decomposed.vcf.gz --threads ${THREADS} -R PRS_vars_chr${CHR}.bed > PRS_chr${CHR}.vcf.gz
plink --vcf PRS_chr${CHR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out PRS_chr${CHR}
done

#####################################
## create plink files for CCSS exp ##
#####################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome
for i in {1..22}; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
unset VCF; \
export THREADS=4; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/"; \
bsub \
        -P "${CHR}_plk" \
        -J "${CHR}_plk" \
        -o "${WORKDIR}/logs/${VCF%.vcf*}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "./plinkQC.sh"; \
done;

#!/usr/bin/bash
module load plink/1.90b
plink --vcf "${WORKDIR}/CCSS_exp_biallelic_${CHR}_ID_edited.vcf.gz" --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out "${WORKDIR}/CCSS_exp_biallelic_${CHR}_ID_edited"

#####################################
## create plink files for CCSS org ##
#####################################
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/
for i in {1..22}; do \
export CHR="chr${i}"; \
echo "splitting $CHR"; \
export THREADS=4; \
export WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/"; \
export VCFDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/"; \
bsub \
        -P "${CHR}_plk" \
        -J "${CHR}_plk" \
        -o "${WORKDIR}/logs/${CHR}_biallelic.%J" \
        -n ${THREADS} \
        -R "rusage[mem=8192]" \
        "./plinkQC.sh"; \
done;


#!/usr/bin/bash
module load plink/1.90b
plink --vcf "${VCFDIR}/${CHR}.dose.vcf.gz" --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out "${WORKDIR}/CCSS_org_GRCh37_${CHR}"



# #!/usr/bin/bash
# module load bcftools/1.9
# module load tabix
# cd ${WORKDIR}
# /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/liftover/liftOver "${VCFDIR}/${CHR}.dose.vcf.gz" /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/tools/liftover/hg19ToHg38.over.chain ${CHR}_GRCh38.dose.vcf.gz ${CHR}_unmapped.bed
# tabix -p vcf "${WORKDIR}/${CHR}_GRCh38.dose.vcf.gz"


## Rename all BIM files with CHR:POS:REF:ALT in shell script in ccss_org
for i in {1..22}; do \
CHR="chr${i}";
awk '{print "chr"$1, "chr"$2":"$6":"$5, $3, $4, $5, $6}' CCSS_org_GRCh37_${CHR}.bim > tmp_bim
mv tmp_bim CCSS_org_GRCh37_${CHR}.bim
done