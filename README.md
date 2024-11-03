# VTA
A Complete Library for Transcriptome Analysis

## Setup

Compile the files:
```sh
make
```

Install Python packages:
```sh
pip install -r requirements.txt
```

## Raw Data Processing

Index the reference genome:
```sh
./index <GRCh38.fna>
```

Align the sequences:
```sh
./align <GRCh38.fna> <input1.fq> <input2.fq> <output.sam>
```

Profile the mRNA expression:
```sh
./profile <gencode.gtf> <output.sam> <output.tsv>
```

*(This library mainly focuses on the analysis of aligned data. If you found the index and align parts too slow, you can use other tools such as [bwa](https://github.com/lh3/bwa) for alignment.)*

## Expression Analysis

Summarize the profiles:
```sh
python3 summarize.py <cohort.csv> <profile_dir/> <start_gene> <value_type>
```

Perform clustering analysis:
```sh
python3 cluster.py <cohort/summary.csv> <group_column1> <group_column2> ... 
```

Perform differentially expressed gene analysis:
```sh
python3 DEG.py <cohort/summarize.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... 
```

Perform gene set enrichment analysis:
```sh
python3 GSEA.py <cohort/summarize.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... <geneset.gmt>
```