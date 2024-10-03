# wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
# tar -zxvf sratoolkit.current-ubuntu64.tar.gz
# cd sratoolkit.3.X.X-ubuntu64/bin/

# Read each line from the file
while read -r gsm_id sra_id; do
    # Check if the output directory already exists
    if [ -d "$sra_id" ]; then
        echo "Folder $sra_id already exists. Skipping..."
        continue
    fi

    # Prefetch the SRA data
    echo "Downloading $sra_id..."
    ./prefetch "$sra_id"

    # Create output directory
    # mkdir -p "$sra_id"

    # Dump the SRA data to the specified directory
    echo "Converting $sra_id to FASTQ..."
    ./fastq-dump --split-files --outdir "$sra_id" "$sra_id/$sra_id.sra"

    echo "$sra_id processing complete."

done < liv_list.txt
