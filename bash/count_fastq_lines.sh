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
echo "${file1} has ${file1_lines} lines; ${file2} has ${file2_lines} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/read_counts_HLHS.txt
done


./GW3-7/V350082588_L03_B5GHUMqzinRABCA-570_1.fq.gz has 87524113 lines; ./GW3-7/V350082588_L03_B5GHUMqzinRABCA-570_2.fq.gz has 167690820 lines
./GW21-6/V350095281_L01_B5GHUMqzinRABEA-582_1.fq.gz has 277092808 lines; ./GW21-6/V350095281_L01_B5GHUMqzinRABEA-582_2.fq.gz has 99973927 lines
./GW22-36/V350095281_L02_B5GHUMqzinRABFA-503_1.fq.gz has 84899558 lines; ./GW22-36/V350095281_L02_B5GHUMqzinRABFA-503_2.fq.gz has 0 lines
./GW22-36/V350095281_L02_B5GHUMqzinRABFA-507_1.fq.gz has 0 lines; ./GW22-36/V350095281_L02_B5GHUMqzinRABFA-507_2.fq.gz has 203424196 lines
./GW22-36/V350095281_L02_B5GHUMqzinRABFA-508_1.fq.gz has 0 lines; ./GW22-36/V350095281_L02_B5GHUMqzinRABFA-508_2.fq.gz has 175942060 lines
./GW124-9/V350095281_L03_B5GHUMqzinRABGA-535_1.fq.gz has 89095293 lines; ./GW124-9/V350095281_L03_B5GHUMqzinRABGA-535_2.fq.gz has 0 lines
./GW132-2B/V350095490_L01_B5GHUMqzinRABHA-542_1.fq.gz has 137866796 lines; ./GW132-2B/V350095490_L01_B5GHUMqzinRABHA-542_2.fq.gz has 0 lines


ls ./*//GW21-6/V350095281_L01_B5GHUMqzinRABEA-582_*.fq.gz | parallel 'echo {} has $(zcat {} | wc -l) lines'

ls ./*/GW3-7/V350082588_L03_B5GHUMqzinRABCA-570_*.fq.gz | parallel 'echo {} has $(zcat {} | wc -l) lines'


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