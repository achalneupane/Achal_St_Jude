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