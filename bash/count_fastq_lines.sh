for file1 in $(ls ./*/*_1.fq.gz| sort -V); do
file2="$(echo $file1 | sed 's/_1.fq.gz/_2.fq.gz/g')"
file1_lines="$(zcat ${file1}| wc -l)"
file2_lines="$(zcat ${file2}| wc -l)"
echo "${file1} has ${file1_lines} lines; ${file2} has ${file2_lines} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/read_counts.txt
done



Re-download:

 # zcat ./*/*/V300063238_L2_B5GHUMrcaRAAADAAA-548_1.fq.gz| wc -l; zcat ./*/*/V300063238_L2_B5GHUMrcaRAAADAAA-548_2.fq.gz| wc -l

./sample/Fastq_Read1; ./sample/Fastq_Read2



for file1 in $(ls ./*/*_1.fq.gz| sort -V); do
file2="$(echo $file1 | sed 's/_1.fq.gz/_2.fq.gz/g')"
file1_lines="$(zcat ${file1}| wc -l)"
file2_lines="$(zcat ${file2}| wc -l)"
echo "${file1} has ${file1_lines} lines; ${file2} has ${file2_lines} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/read_counts_Ashley.txt
done





get_counts () {
file1=$1
file2="$(echo $file1 | sed 's/_1.fq.gz/_2.fq.gz/g')"
# file1_lines="$(zcat ${file1}| wc -l)"
# file2_lines="$(zcat ${file2}| wc -l)"
# echo "${file1} has ${file1_lines} lines; ${file2} has ${file2_lines} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/all_read_counts_from_function.txt
echo "${file1} lines; ${file2} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/all_read_counts_from_function.txt
}

export -f get_counts


parallel -j8 get_counts {} ::: $(ls ./*/*/*_1.fq.gz| sort -V)

parallel dowork ::: "${list[@]}" ::: "${other[@]}"

ls ./*/*/*_1.fq.gz| sort -V



ls ./*/*/*_1.fq.gz| sort -V | parallel -j30 'var=$(echo {}|sed s/_1.fq.gz/_2.fq.gz/g); echo {} has $(zcat -n {} | wc -l) lines and $var has $(zcat $var | wc -l)'  >> all_fastq_line_counts_v2.txt

for file in $(ls *.fq.gz| sort -V); do
md5sum $file
done


# gzip: V300009670_L4_B5GHUMsgfRAAALAAA-522_1.fq.gz: not in gzip format
# gzip: V300009670_L4_B5GHUMsgfRAAALAAA-522_2.fq.gz: not in gzip format



for dir in */ ; do     echo $dir | wc -l; done


V300009599_L3_B5GHUMsgfRAAAZAAA-503_1.fq.gz



1607 chr10

994 chr2


## /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_1kgp/look_up_unmatched_variants.sh
## Check variants in 1000G file: /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_1kgp
cat breast_cancer_prs_vars_Mavaddat_2019_ER_OVERALL_Breast_not_matched.txt| cut -d$'\t' -f2| parallel -j30 'chr=$(grep {} breast_cancer_prs_vars_Mavaddat_2019_ER_OVERALL_Breast_not_matched.txt|sed "s/\s.*$//"); zcat CCSS.AnalysisSet.${chr}.gen.gz| grep -w {}' >> matched_positions.txt
cat /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/Mavaddat_2019_ER_OVERALL_Breast_not_found.txt | while read line; do
echo ${line}
grep -w ${line} matched_positions.txt| wc -l
done
## Zero match from the ones missing


cat breast_cancer_prs_vars_Mavaddat_2019_ER_OVERALL_Breast_not_matched.txt | while read line; do
var="$(echo $line| cut -d' ' -f2)"
echo $var
# done
grep -w ${var} matched_positions.txt| wc -l
done

73 did not match in Breast cancer


# Also check for sarcoma, meningioma, pleitropy. Since fewer breast cancer variants were found compared to ccss_original gwas data, I am only searching for the ones not found here
sarcoma_meningioma_pleiotropy_not_matched.txt
22      44924073        A       G       0.0134  Mavaddat_2019_ER_OVERALL_Breast Breast  Y       chr22:44924073  45319953        22:45319953     22:45319943-45319963
20      63636980        G       A       0.002199633     Pleiotropy_PRSWEB       Pleiotropy      Y       chr20:63636980  62268333        20:62268333     20:62268323-62268343
20      63636988        T       C       0.000191582     Pleiotropy_PRSWEB       Pleiotropy      Y       chr20:63636988  62268341        20:62268341     20:62268331-62268351
20      63640895        C       T       0.0012619       Pleiotropy_PRSWEB       Pleiotropy      Y       chr20:63640895  62272248        20:62272248     20:62272238-62272258
20      63686867        T       C       0.001131128     Pleiotropy_PRSWEB       Pleiotropy      Y       chr20:63686867  62318220        20:62318220     20:62318210-62318230
20      63695226        G       T       0.000175282     Pleiotropy_PRSWEB       Pleiotropy      Y       chr20:63695226  62326579        20:62326579     20:62326569-62326589
20      63711343        T       C       0.000171664     Pleiotropy_PRSWEB       Pleiotropy      Y       chr20:63711343  62342695        20:62342695     20:62342685-62342705
20      64071667        T       C       0.000278223     Pleiotropy_PRSWEB       Pleiotropy      Y       chr20:64071667  62703020        20:62703020     20:62703010-62703030
11      258909  A       T       0.364643114     Meningioma_Claus        Meningioma      Y       chr11:258909    258909  11:258909       11:258899-258919
1       10986508        C       T       0.717839793     Sarcoma_Machiela        Sarcoma Y       chr1:10986508   11046565        1:11046565      1:11046555-11046575


cat sarcoma_meningioma_pleiotropy_not_matched.txt| cut -d$'\t' -f2| parallel -j30 'chr=$(grep {} sarcoma_meningioma_pleiotropy_not_matched.txt|sed "s/\s.*$//"); zcat CCSS.AnalysisSet.${chr}.gen.gz| grep -w {}' >> matched_positions_sarcoma_meningioma_pleitropy_not_matched.txt
