#!/bin/bash
module load pigz
module load conda3/202210

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy_round2/newvolume/expdata/Sub_library_${NUM}

# zcat Sub_library_${NUM}_*_1.fq.gz| pigz -p ${THREADS} > S${NUM}_R1.fq.gz
# zcat Sub_library_${NUM}_*_*_2.fq.gz| pigz -p ${THREADS} > S${NUM}_R2.fq.gz

conda activate /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/envs/spipe
split-pipe \
    --mode all \
    --chemistry v2 \
    --nthreads ${THREADS} \
    --genome_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/newvolume/genomes/hg38/ \
    --fq1 /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy_round2/newvolume/expdata/Sub_library_${NUM}/S${NUM}_R1.fq.gz \
    --output_dir /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy_round2/newvolume/analysis/S${NUM}-out \
    --samp_sltab /research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy_round2/newvolume/analysis/WT_100K_Sample_Loading_Table.xlsm

# End of script

## To submit
# for iter in {1..8}; do \
#   unset NUM; \
#   export NUM=${iter}; \
#   echo "Doing Sublibrary ${NUM}"; \
#   export MEM=12; \
#   export THREADS=12; \
#   bsub \
#   -P "Sublib${NUM}" \
#   -J "Sublib${NUM}" \
#   -o "Sublib${NUM}_e.%J" \
#   -n ${THREADS} \
#   -R "rusage[mem=8192]" \
#   "/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy_round2/newvolume/expdata/sub.bsub.sh"; \
# done; 