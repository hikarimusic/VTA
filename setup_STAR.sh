#!/bin/bash
set -e
THREADS=8
WORK_DIR="$PWD"

# Function to download
download_file() {
    local url="$1"
    local output="$2"
    if ! wget -q --show-progress -O "$output" "$url"; then
        echo "[Error] Cannot download $url"
        return 1
    fi
}

# 1. Download STAR
if [ ! -f "${WORK_DIR}/STAR" ]; then
    echo "[STAR] Downloading ..."
    download_file "https://raw.githubusercontent.com/alexdobin/STAR/master/bin/Linux_x86_64_static/STAR" "${WORK_DIR}/STAR"
    chmod +x "${WORK_DIR}/STAR"
else
    echo "[STAR] Complete"
fi

# 2. Download genome
GENOME_FILE="${WORK_DIR}/GRCh38.primary_assembly.genome.fa"
if [ ! -f "${GENOME_FILE}" ]; then
    echo "[Genome] Downloading ..."
    download_file "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz" "${GENOME_FILE}.gz"
    gunzip "${GENOME_FILE}.gz"
else
    echo "[Genome] Complete"
fi

# 3. Download gencode
GTF_FILE="${WORK_DIR}/gencode.v36.annotation.gtf"
if [ ! -f "${GTF_FILE}" ]; then
    echo "[Gencode] Downloading ..."
    download_file "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz" "${GTF_FILE}.gz"
    gunzip "${GTF_FILE}.gz"
else
    echo "[Gencode] Complete"
fi

# 4. Build index
STAR_INDEX_DIR="${WORK_DIR}/STAR_index"
if [ ! -f "${STAR_INDEX_DIR}/SAindex" ]; then
    echo "[Index] Building ..."
    mkdir -p ${STAR_INDEX_DIR}
    AVAILABLE_RAM=$(free -b | awk '/Mem:/ {print $2}')
    LIMIT_RAM=$(echo "${AVAILABLE_RAM} * 0.9" | bc | cut -d'.' -f1)
    ./STAR --runMode genomeGenerate \
           --runThreadN ${THREADS} \
           --genomeDir ${STAR_INDEX_DIR} \
           --genomeFastaFiles ${GENOME_FILE} \
           --sjdbGTFfile ${GTF_FILE} \
           --sjdbOverhang 99 \
           --limitGenomeGenerateRAM ${LIMIT_RAM}
else
    echo "[Index] Complete"
fi
