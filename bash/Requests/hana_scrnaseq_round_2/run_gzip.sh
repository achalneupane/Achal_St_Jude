base_dir="/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy_round2/newvolume/expdata"
# Function to process a sub-library
process_sub_library() {
    sub_dir="$1"
    cd "${sub_dir}"

    # Print a message indicating the start of processing
    echo "Processing sub-library in: $sub_dir"

    zcat Sub_library_*_1.fq.gz | gzip -c > "S${i}_R1.fq.gz"
    zcat Sub_library_*_2.fq.gz | gzip -c > "S${i}_R2.fq.gz"

    # Print a message indicating the completion of processing
    echo "Finished processing sub-library in: $sub_dir"
}

# Number of parallel processes (adjust as needed)
num_processes=8
# Array to store background process IDs
process_ids=()
# Loop through the sub-libraries and process them in parallel
for i in {1..8}
do
    sub_dir="${base_dir}/Sub_library_${i}"
    process_sub_library "${sub_dir}" &
    process_ids+=($!)

    # Limit the number of parallel processes
    if [[ ${#process_ids[@]} -ge $num_processes ]]; then
        for pid in "${process_ids[@]}"; do
            wait "$pid"
        done
        process_ids=()
    fi
done
# Wait for all background processes to finish
for pid in "${process_ids[@]}"; do
    wait "$pid"
done
