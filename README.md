# VTA
A Complete Library for NGS Data Analysis

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
python3 summarize.py <cohort_file.csv> <profile_dir/> <start_gene> <value_type>
```

Generate PCA plots:
```sh
python3 PCA.py <cohort_file.csv> <group_column>
```

Perform clustering analysis:
```sh
python3 clustering.py <cohort_file.csv> <group_column>
```

Analyze differentially expressed gene:
```sh
python3 DEG.py <cohort_file.csv> <group_column> <group1-a> <group1-b> ... -- <group2-a> <group2-b> ...
```