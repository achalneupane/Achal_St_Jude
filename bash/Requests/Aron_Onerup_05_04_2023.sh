# Request from Aron Onerup emailed on 05/04/2023
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned
cd Aron_Onerup_05_04_2023

RSID
rs1056892
rs1786814
rs2232228
rs2229774
rs17863783
rs7853758
rs4982753
rs4149178

##################
## ccss_exp_WGS ##
##################

# GrCh38
cat <<\EoF > extract_vars_GRCH38.bed
#chrom chromStart chromEnd name ref alts
chr5 609977 609978 rs924607 C T,
EoF


plink --bfile ../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf.gz --chr 5 --from-bp 609977 --to-bp 609978 --make-bed --out aron_bfile
# plink --bfile ../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.preQC_biallelic_renamed_ID_edited.vcf.gz --snp chr5:609978:C:T --make-bed --out mydata_chr5:609978:C:T

plink --bfile aron_bfile --recodeA --out aron_bfile_recodeA 

zcat /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr5.PASS.decomposed.vcf.gz |head -5000|  grep "#CHROM" | tr "\t" "\n " | tail -n +10 | uniq > QCed_samples

