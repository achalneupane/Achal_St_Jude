# Request from Yutaka emailed on 09/27/2022


RSID
rs1056892
rs1786814
rs2232228
rs2229774
rs17863783
rs7853758
rs4982753
rs4149178

# GrCh38
cat <<\EoF > extract_vars_GRCH38.bed
#chrom chromStart chromEnd name ref alts
chr2 233693630 233693631 rs17863783 G T,
chr6 43304449 43304450 rs4149178 A G,
chr9 84286010 84286011 rs7853758 G A,
chr12 53211760 53211761 rs2229774 G A,
chr14 23345359 23345360 rs4982753 C A,T,
chr16 69109673 69109674 rs2232228 A C,G,
chr18 37497064 37497065 rs1786814 G A,C,
chr21 36146407 36146408 rs1056892 G A,
EoF

# GRCh37
cat <<\EoF > extract_vars_GRCH37.bed
#chrom chromStart chromEnd name ref alts
chr2 234602276 234602277 rs17863783 G T,
chr6 43272187 43272188 rs4149178 A G,
chr9 86900925 86900926 rs7853758 G A,
chr12 53605544 53605545 rs2229774 G A,
chr14 23814568 23814569 rs4982753 C A,T,
chr16 69143576 69143577 rs2232228 A C,G,
chr18 35077027 35077028 rs1786814 G A,C,
chr21 37518705 37518706 rs1056892 G A,
EoF

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022
ln -s ../CCSS_exp_biallelic*vcf.gz* .

module load bcftools/1.9
module load plink/1.90b

cat extract_vars_GRCH38.txt |sed 's/^[ \t]*//;s/[ \t]*$//' > extract_vars_GRCH38_edited.txt

IFS=$'\n' # set IFS
for line in $(cat extract_vars_GRCH38.bed| sed '1d'); do
VAR="$(echo ${line}| awk '{print $1":"$2"-"$3}')"
CHR="$(echo $VAR |awk -F':' '{print $1}')"
echo "Doing ${VAR}"
bcftools view /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_${CHR}_ID_edited.vcf.gz ${VAR} > /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data/yutaka_${VAR}.vcf.gz
bcftools query -f '%CHROM:%POS:%REF:%ALT[\n%SAMPLE=%GT]\n' /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data/yutaka_${VAR}.vcf.gz >> /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data/yutaka_${VAR}.geno.txt
plink --vcf /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data/yutaka_${VAR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --make-bed --out /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data/yutaka_${VAR} 2>&1 | tee -a extract_plink_all.log
done

# paste *.geno.txt | sed 's/^[ \t]*//;s/[ \t]*$//' > all_genotypes_from_ccss_WGS_for_Yutaka.txt

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data
ls *.bim| sort -V | sed 's/\.bim//g'|sed -n '1d;p' > merge_list.list

plink --bfile yutaka_chr2:233693630-233693631 --merge-list merge_list.list --keep-allele-order --out merged.dat.yutaka

plink --bfile merged.dat.yutaka --recodeA --out merged.dat.yutaka_recodeA 

