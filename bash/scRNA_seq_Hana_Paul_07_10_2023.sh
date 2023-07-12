cd /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_1
## Sublibrary_1; concat R1 and R2 in different lanes
zcat Sub_library_1_CKDL230017038-1A*_L*_1* | gzip -c > S1_R1.fq.gz
zcat Sub_library_1_CKDL230017038-1A*_L*_2* | gzip -c > S1_R2.fq.gz

split-pipe \
    --mode all \
    --chemistry v2 \
    --nthreads 12 \
    --genome_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/genomes/hg38/ \
    --fq1 /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_1/S1_R1.fq.gz \
    --output_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S1-out \
    --samp_sltab /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/Hana_WT_mini.xlsm


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_2
## Sublibrary_1; concat R1 and R2 in different lanes
zcat Sub_library_2_CKDL230017039-1A_*_L*_1.fq.gz | gzip -c > S2_R1.fq.gz
zcat Sub_library_2_CKDL230017039-1A_*_L*_2.fq.gz | gzip -c > S2_R2.fq.gz

split-pipe \
    --mode all \
    --chemistry v2 \
    --nthreads 12 \
    --genome_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/genomes/hg38/ \
    --fq1 /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/expdata/Sub_library_2/S2_R1.fq.gz \
    --output_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S2-out \
    --samp_sltab /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/Hana_WT_mini.xlsm



    
## Combine sublibraries
split-pipe \
    --mode comb \
    --sublibraries /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S1-out \
                   /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S2-out \
    --output_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/analysis/S_combined




