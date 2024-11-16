#!/bin/bash
set -e

cat > process_star_output.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
from collections import defaultdict

def parse_gtf(gtf_file):
    gene_exons = defaultdict(lambda: [])
    gene_info = {}
    
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            
            # Parse attributes
            attr_string = fields[8]
            attrs = {}
            for attr in attr_string.strip().split(';'):
                if attr.strip():
                    try:
                        key, value = attr.strip().split(' ', 1)
                        attrs[key] = value.strip('"')
                    except ValueError:
                        continue
            if 'gene_id' not in attrs:
                continue
            gene_id = attrs['gene_id']
            
            # Store gene information
            if fields[2] == 'gene':
                if 'gene_name' in attrs and 'gene_type' in attrs:
                    gene_info[gene_id] = {
                        'gene_name': attrs['gene_name'],
                        'gene_type': attrs['gene_type']
                    }
            
            # Store exon coordinates
            elif fields[2] == 'exon':
                start = int(fields[3])
                end = int(fields[4])
                gene_exons[gene_id].append((start, end))
    
    # Calculate gene length
    gene_lengths = {}
    for gene_id, exons in gene_exons.items():
        if not exons:
            continue
        
        # Merge overlapping exons
        sorted_exons = sorted(exons)
        union_length = 0
        current_start, current_end = sorted_exons[0]

        for start, end in sorted_exons[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                union_length += current_end - current_start + 1
                current_start, current_end = start, end
        
        union_length += current_end - current_start + 1
        gene_lengths[gene_id] = union_length
    
    return gene_lengths, gene_info

def calculate_expression_metrics(counts_file, gene_lengths, gene_info, output_file):
    # Gene count
    df = pd.read_csv(counts_file, sep='\t', header=None, skiprows=4,
                     names=['gene_id', 'unstranded', 'stranded_first', 'stranded_second'])
    
    df['length'] = df['gene_id'].map(gene_lengths)
    df['gene_name'] = df['gene_id'].map(lambda x: gene_info.get(x, {}).get('gene_name', ''))
    df['gene_type'] = df['gene_id'].map(lambda x: gene_info.get(x, {}).get('gene_type', ''))
    
    # Calculate TPM
    df['rpk_unstranded'] = df['unstranded'] / (df['length'] / 1000)
    total_rpk = df['rpk_unstranded'].sum()
    df['tpm_unstranded'] = df['rpk_unstranded'] / total_rpk * 1e6
    
    # Calculate FPKM
    total_reads = df['unstranded'].sum()
    df['fpkm_unstranded'] = df['unstranded'] * 1e9 / (df['length'] * total_reads)
    
    # Calculate FPKM-UQ
    uq = df[df['unstranded'] > 0]['unstranded'].quantile(0.75)
    df['fpkm_uq_unstranded'] = df['unstranded'] * 1e9 / (df['length'] * uq * len(gene_lengths))
    
    # Header
    header_stats = pd.read_csv(counts_file, sep='\t', header=None, nrows=4,
                              names=['gene_id', 'unstranded', 'stranded_first', 'stranded_second'])
    header_stats_df = pd.DataFrame([
        ['N_unmapped', '', '', header_stats.iloc[0, 1], header_stats.iloc[0, 2], header_stats.iloc[0, 3], '', '', ''],
        ['N_multimapping', '', '', header_stats.iloc[1, 1], header_stats.iloc[1, 2], header_stats.iloc[1, 3], '', '', ''],
        ['N_noFeature', '', '', header_stats.iloc[2, 1], header_stats.iloc[2, 2], header_stats.iloc[2, 3], '', '', ''],
        ['N_ambiguous', '', '', header_stats.iloc[3, 1], header_stats.iloc[3, 2], header_stats.iloc[3, 3], '', '', '']
    ], columns=['gene_id', 'gene_name', 'gene_type', 'unstranded', 'stranded_first', 'stranded_second', 
                'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded'])
    
    # Final output
    result_df = df[['gene_id', 'gene_name', 'gene_type', 'unstranded', 'stranded_first', 'stranded_second',
                    'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']]
    result_df = result_df.sort_values('gene_id')
    final_df = pd.concat([header_stats_df, result_df], ignore_index=True)
    with open(output_file, 'w') as f:
        f.write('# gene-model: GENCODE v36\n')
        final_df.to_csv(f, sep='\t', index=False)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: script.py <ReadsPerGene.out.tab> <gencode.gtf> <output.tsv>")
        sys.exit(1)
        
    counts_file = sys.argv[1]
    gtf_file = sys.argv[2]
    output_file = sys.argv[3]
    
    gene_lengths, gene_info = parse_gtf(gtf_file)
    
    calculate_expression_metrics(counts_file, gene_lengths, gene_info, output_file)
EOF

INPUT_DIR="$1"
OUTPUT_DIR="$2"
WORK_DIR="$PWD"
THREADS=16

mkdir -p "${OUTPUT_DIR}"

# Run STAR and process output
process_reads() {
    local base_id="$1"
    local reads="$2"
    local output_prefix="${OUTPUT_DIR}/${base_id}"
    
    echo "[Mapping] ${base_id}"
    
    # Run STAR alignment
    ./STAR --runThreadN ${THREADS} \
           --genomeDir "${WORK_DIR}/STAR_index" \
           --readFilesIn ${reads} \
           --outFileNamePrefix "${output_prefix}." \
           --quantMode GeneCounts
    
    echo "[Transform] ${base_id}"
    
    python3 "${WORK_DIR}/process_star_output.py" \
            "${output_prefix}.ReadsPerGene.out.tab" \
            "${WORK_DIR}/gencode.v36.annotation.gtf" \
            "${output_prefix}.tsv"

    # Clean up intermediate files
    for file in "${output_prefix}."*; do
        if [ "${file}" != "${output_prefix}.tsv" ]; then
            rm -f "${file}"
        fi
    done
}

# Process FASTQ files
find "${INPUT_DIR}" -name "*.fastq" | while read -r fastq; do
    base_path="${fastq%_[12].fastq}"
    base_id=$(basename "${base_path}")
    
    # paired-end
    if [[ "${fastq}" =~ _1\.fastq$ ]]; then
        fastq2="${base_path}_2.fastq"
        if [ -f "${fastq2}" ]; then
            process_reads "${base_id}" "${fastq} ${fastq2}"
            continue
        fi
    fi
    
    # single-end
    if [[ ! "${fastq}" =~ _2\.fastq$ ]]; then
        process_reads "${base_id}" "${fastq}"
    fi
done

rm -f "${WORK_DIR}/process_star_output.py"