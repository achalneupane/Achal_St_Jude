for file1 in $(ls ./*/*_1.fq.gz| sort -V); do
file2="$(echo $file1 | sed 's/_1.fq.gz/_2.fq.gz/g')"
file1_lines="$(zcat ${file1}| wc -l)"
file2_lines="$(zcat ${file2}| wc -l)"
echo "${file1} has ${file1_lines} lines; ${file2} has ${file2_lines} lines " >> /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/read_counts.txt
done



Re-download:

V300009702_L1_B5GHUMsgfRAAAGAAA-511_1.fq.gz of JR44-C4-p21 = Not OK
V100004093_L1_B5GHUMqluRAABHAAA-529_1.fq.gz of SC08-C1-p18 = Not OK
V100004093_L2_B5GHUMqluRAABHAAA-532_1.fq.gz of SC08-C1-p18 = Not OK
V100004093_L3_B5GHUMqluRAABIAAA-539_1.fq.gz of SC10-C14-p12 = Not OK
V300009599_L3_B5GHUMsgfRAAAZAAA-503_1.fq.gz of SC26-C8-p9 = Not OK
V300009599_L3_B5GHUMsgfRAAAZAAA-508_1.fq.gz of SC26-C8-p9 = Not OK

V300063449_L3_B5GHUMrcaRAAAAAAA-515_1.fq.gz of JW1-12 = Not OK

 V300063238_L2_B5GHUMrcaRAAADAAA-548_1.fq.gz



./JW4-5/V300063238_L2_B5GHUMrcaRAAADAAA-548_1.fq.gz has 86235133 lines; ./JW4-5/V300063238_L2_B5GHUMrcaRAAADAAA-548_2.fq.gz has 99726544 lines
./JR44-C4-p21/V300009702_L1_B5GHUMsgfRAAAGAAA-511_1.fq.gz has 108811744 lines; ./JR44-C4-p21/V300009702_L1_B5GHUMsgfRAAAGAAA-511_2.fq.gz has 108712383 lines
./SC08-C1-p18/V100004093_L1_B5GHUMqluRAABHAAA-529_1.fq.gz has 50991480 lines; ./SC08-C1-p18/V100004093_L1_B5GHUMqluRAABHAAA-529_2.fq.gz has 149033756 lines
./SC08-C1-p18/V100004093_L2_B5GHUMqluRAABHAAA-532_1.fq.gz has 71879828 lines; ./SC08-C1-p18/V100004093_L2_B5GHUMqluRAABHAAA-532_2.fq.gz has 36752739 lines
./SC10-C14-p12/V100004093_L3_B5GHUMqluRAABIAAA-539_1.fq.gz has 21693855 lines; ./SC10-C14-p12/V100004093_L3_B5GHUMqluRAABIAAA-539_2.fq.gz has 125494356 lines
./SC26-C8-p9/V300009599_L3_B5GHUMsgfRAAAZAAA-503_1.fq.gz has 20699537 lines; ./SC26-C8-p9/V300009599_L3_B5GHUMsgfRAAAZAAA-503_2.fq.gz has 117567984 lines
./SC26-C8-p9/V300009599_L3_B5GHUMsgfRAAAZAAA-508_1.fq.gz has 19952443 lines; ./SC26-C8-p9/V300009599_L3_B5GHUMsgfRAAAZAAA-508_2.fq.gz has 163416120 lines
