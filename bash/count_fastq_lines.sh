for file1 in $(ls ./*/*_1.fq.gz| sort -V); do
file2="$(echo $file1 | sed 's/_1.fq.gz/_2.fq.gz/g')"
file1_lines="$(zcat ${file1}| wc -l)"
file2_lines="$(zcat ${file2}| wc -l)"
echo "${file1} has ${file1_lines} lines; ${file2} has ${file2_lines} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/read_counts.txt
done



Re-download:

 # zcat ./*/*/V300063238_L2_B5GHUMrcaRAAADAAA-548_1.fq.gz| wc -l; zcat ./*/*/V300063238_L2_B5GHUMrcaRAAADAAA-548_2.fq.gz| wc -l

./sample/Fastq_Read1; ./sample/Fastq_Read2

LMNA
./JW4-5/V300063238_L2_B5GHUMrcaRAAADAAA-548_1.fq.gz has 86235133 lines; ./JW4-5/V300063238_L2_B5GHUMrcaRAAADAAA-548_2.fq.gz has 99726544 lines


LVNC
./GW115-8/V350095253_L02_B5GHUMqzinRAAUB-569_1.fq.gz has 142948056 lines; ./GW115-8/V350095253_L02_B5GHUMqzinRAAUB-569_2.fq.gz has 160273864 lines
./GW129-24/V350097135_L01_B5GHUMobpeRAAIA-582_1.fq.gz has 133064759 lines; ./GW129-24/V350097135_L01_B5GHUMobpeRAAIA-582_2.fq.gz has 0 lines
./GW159-11/V350095253_L03_B5GHUMqzinRAAWB-581_1.fq.gz has 243669836 lines; ./GW159-11/V350095253_L03_B5GHUMqzinRAAWB-581_2.fq.gz has 185581124 lines
./GW159-11/V350095253_L03_B5GHUMqzinRAAWB-582_1.fq.gz has 137322658 lines; ./GW159-11/V350095253_L03_B5GHUMqzinRAAWB-582_2.fq.gz has 276421320 lines
./GW159-11/V350095253_L03_B5GHUMqzinRAAWB-585_1.fq.gz has 112984799 lines; ./GW159-11/V350095253_L03_B5GHUMqzinRAAWB-585_2.fq.gz has 66228385 lines
./GW165-11/V350095253_L04_B5GHUMqzinRAAXA-503_1.fq.gz has 196076684 lines; ./GW165-11/V350095253_L04_B5GHUMqzinRAAXA-503_2.fq.gz has 0 lines
./GW167-6/V350095186_L04_B5GHUMqzinRDAAA-551_1.fq.gz has 0 lines; ./GW167-6/V350095186_L04_B5GHUMqzinRDAAA-551_2.fq.gz has 200450684 lines


Previous download
./JR44-C4-p21/V300009702_L1_B5GHUMsgfRAAAGAAA-511_1.fq.gz has 108811744 lines; ./JR44-C4-p21/V300009702_L1_B5GHUMsgfRAAAGAAA-511_2.fq.gz has 108712383 lines
./SC08-C1-p18/V100004093_L1_B5GHUMqluRAABHAAA-529_1.fq.gz has 50991480 lines; ./SC08-C1-p18/V100004093_L1_B5GHUMqluRAABHAAA-529_2.fq.gz has 149033756 lines
./SC08-C1-p18/V100004093_L2_B5GHUMqluRAABHAAA-532_1.fq.gz has 71879828 lines; ./SC08-C1-p18/V100004093_L2_B5GHUMqluRAABHAAA-532_2.fq.gz has 36752739 lines
./SC10-C14-p12/V100004093_L3_B5GHUMqluRAABIAAA-539_1.fq.gz has 21693855 lines; ./SC10-C14-p12/V100004093_L3_B5GHUMqluRAABIAAA-539_2.fq.gz has 125494356 lines
./SC26-C8-p9/V300009599_L3_B5GHUMsgfRAAAZAAA-503_1.fq.gz has 20699537 lines; ./SC26-C8-p9/V300009599_L3_B5GHUMsgfRAAAZAAA-503_2.fq.gz has 117567984 lines
./SC26-C8-p9/V300009599_L3_B5GHUMsgfRAAAZAAA-508_1.fq.gz has 19952443 lines; ./SC26-C8-p9/V300009599_L3_B5GHUMsgfRAAAZAAA-508_2.fq.gz has 163416120 lines









zcat ./*/*/V350095508_L03_B5GHUMqzinRAAOA-582_1.fq.gz| wc -l; zcat ./*/*/V350095508_L03_B5GHUMqzinRAAOA-582_2.fq.gz| wc -l
zcat ./*/*/V350095508_L03_B5GHUMqzinRAAOA-585_1.fq.gz| wc -l; zcat ./*/*/V350095508_L03_B5GHUMqzinRAAOA-585_2.fq.gz| wc -l


grep V350095186_L04_B5GHUMqzinRDAAA-551_1.fq.g ../../*.txt


for file1 in $(ls ./*/*_1.fq.gz| sort -V); do
file2="$(echo $file1 | sed 's/_1.fq.gz/_2.fq.gz/g')"
file1_lines="$(zcat ${file1}| wc -l)"
file2_lines="$(zcat ${file2}| wc -l)"
echo "${file1} has ${file1_lines} lines; ${file2} has ${file2_lines} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/read_counts_HLHS.txt
done
