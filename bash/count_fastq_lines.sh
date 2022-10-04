for file1 in $(ls ./*/*_1.fq.gz| sort -V); do
file2="$(echo $file1 | sed 's/_1.fq.gz/_2.fq.gz/g')"
file1_lines="$(zcat ${file1}| wc -l)"
file2_lines="$(zcat ${file2}| wc -l)"
echo "${file1} has ${file1_lines} lines; ${file2} has ${file2_lines} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/read_counts.txt
done



Re-download:

V300063238_L2_B5GHUMrcaRAAADAAA-548_1.fq.gz


for file1 in $(ls ./*/*_1.fq.gz| sort -V); do
file2="$(echo $file1 | sed 's/_1.fq.gz/_2.fq.gz/g')"
echo ${file2}
done