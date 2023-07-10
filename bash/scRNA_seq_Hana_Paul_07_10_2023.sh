split-pipe \
    --mode all \
    --chemistry v2 \
    --nthreads 12 \
    --genome_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/genomes/hg38/ \
    --fq1 /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_1/Sub_library_1_CKDL230017038-1A_H73NCDSX7_L1_1.fq.gz \
    --output_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S1-out \
    --samp_sltab /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/Hana_WT_mini.xlsm



split-pipe \
    --mode all \
    --chemistry v2 \
    --nthreads 12 \
    --genome_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/genomes/hg38/ \
    --fq1 /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_1/S1_R1.fq.gz \
    --output_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S1-out \
    --samp_sltab /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/Hana_WT_mini.xlsm




Hi Jason, 

I have a few questions and encountered some issues while running the split-pipe on Sub_library_1.


-rwxrwx--- 1 aneupane Domain Users 8.6G Jul 10 09:46 Sub_library_1_CKDL230017038-1A_H75JGDSX7_L3_2.fq.gz
-rwxrwx--- 1 aneupane Domain Users 7.6G Jul 10 09:46 Sub_library_1_CKDL230017038-1A_H73NCDSX7_L1_2.fq.gz
-rwxrwx--- 1 aneupane Domain Users 9.3G Jul 10 09:46 Sub_library_1_CKDL230017038-1A_H75JGDSX7_L3_1.fq.gz
-rwxrwx--- 1 aneupane Domain Users 8.4G Jul 10 09:47 Sub_library_1_CKDL230017038-1A_H73NCDSX7_L1_1.fq.gz    


Initially, I attempted to concatenate the R1 files from L1 and L3 together, but I was uncertain if this was the correct approach since they come from different runs. Then, I used the following command:

split-pipe \
    --mode all \
    --chemistry v2 \
    --nthreads 12 \
    --genome_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/genomes/hg38/ \
    --fq1 /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_1/S1_R1.fq.gz \
    --output_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S1-out \
    --samp_sltab /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/Hana_WT_mini.xlsm


Additionally, I tried running the command without concatenating the files and used the following command:

split-pipe \
    --mode all \
    --chemistry v2 \
    --nthreads 12 \
    --genome_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/genomes/hg38/ \
    --fq1 /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_1/Sub_library_1_CKDL230017038-1A_H73NCDSX7_L1_1.fq.gz \
    --output_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S1-out \
    --samp_sltab /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/Hana_WT_mini.xlsm


However, in both cases, I encountered the same errors mentioned below. I would appreciate it if you could help me understand where I might be going wrong.

# Init SplitPipe split-pipe v1.0.6p
# Initializing defaults
# Setting threads to 12 (ncpu=64, initial nthreads=12.0) (init threads)
# No random seed given; Creating one
# No fq2; Setting from fq1: /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_1/S1_R2.fq.gz (init fastqs)
# Parsing SampleLoadingTable /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/Hana_WT_mini.xlsm
# --------------------------------------------
# Samples from SampleLoadingTable: /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/Hana_WT_mini.xlsm
#  1 25 A1-A3
#  2 GW64 A4-A6
#  3 GW129 A7-A9
#  4 GW168 A10-A12
# --------------------------------------------
# Table has 4 samples in 12 wells:
#         1  2  3  4  5  6  7  8  9  10  11  12
#  Wells
#  A      1  1  1  2  2  2  3  3  3   4   4   4
# --------------------------------------------
Traceback (most recent call last):
  File "/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/envs/spipe/bin/split-pipe", line 395, in <module>
    ok = main()
  File "/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/envs/spipe/bin/split-pipe", line 190, in main
    spipe = spclass.SplitPipe(VERSION, args)
  File "/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/envs/spipe/lib/python3.9/site-packages/splitpipe/spclass.py", line 126, in __init__
    eval(f"self.{call}")()
  File "/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/envs/spipe/lib/python3.9/site-packages/splitpipe/spclass.py", line 766, in _init_sample_list
    new_name, story = spsamp.clean_sample_def_name(name)
  File "/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/envs/spipe/lib/python3.9/site-packages/splitpipe/spsamp.py", line 686, in clean_sample_def_name
    new_name = utils.clean_vname_chars(name, to_under=True, dash_ok=True)
  File "/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/envs/spipe/lib/python3.9/site-packages/splitpipe/utils.py", line 1140, in clean_vname_chars
    raise Exception(story)
Exception: name not string or list; type: 0
# Closing down SplitPipe
# Total time 0:00:00.42
