# Read each line from the list
while read -r gsm_id sra_id; do
    # Check if the folder exists
    if [ ! -d "/mnt/d/$sra_id" ]; then
        echo "Folder $sra_id does not exist. Skipping..."
        continue
    fi

    # Define the expected output file paths
    sam_file="/mnt/d/$sra_id/$sra_id.sam"
    tsv_file="/mnt/d/$sra_id/$sra_id.tsv"

    # Skip if both .sam and .tsv files already exist
    if [ -f "$sam_file" ] && [ -f "$tsv_file" ]; then
        echo "$sra_id already has .sam and .tsv files. Skipping..."
        continue
    fi

    # Check the existence of the FASTQ files
    fastq_1="/mnt/d/$sra_id/${sra_id}_1.fastq"
    fastq_2="/mnt/d/$sra_id/${sra_id}_2.fastq"
    
    if [ ! -f "$fastq_1" ]; then
        echo "FASTQ file $fastq_1 does not exist. Skipping..."
        continue
    fi

    # Run BWA if the .sam file does not exist
    if [ ! -f "$sam_file" ]; then
        if [ -f "$fastq_2" ]; then
            echo "Running BWA MEM for paired-end data for $sra_id..."
            ./bwa mem GRCh38.NC.fa -t 16 "$fastq_1" "$fastq_2" > "$sam_file"
        else
            echo "Running BWA MEM for single-end data for $sra_id..."
            ./bwa mem GRCh38.NC.fa -t 16 "$fastq_1" > "$sam_file"
        fi
    else
        echo "$sam_file already exists. Skipping BWA MEM..."
    fi

    # Run the profiling command if the .tsv file does not exist
    if [ ! -f "$tsv_file" ]; then
        echo "Profiling $sra_id..."
        ./profile gencode.v46.basic.annotation.gtf "$sam_file" "$tsv_file"
    else
        echo "$tsv_file already exists. Skipping profiling..."
    fi

    echo "$sra_id processing complete."

done < hcc_list.txt
