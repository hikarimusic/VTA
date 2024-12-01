# AMATERASU
Automated Mapping And Transcriptome Expression Research Analysis Suite

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

---
If you want to use STAR to process the data:
```sh
chmod +x setup_STAR.sh process_STAR.sh
```
Download and index genome files:
```sh
./setup_STAR.sh
```
Align and profile RNA-seq data:
```sh
./process_STAR.sh <input_dir/> <output_dir/>
```


## Expression Analysis

Summarize the profiles:
```sh
python3 summarize.py <cohort.csv> <profile_dir/> <value_type>
```

Perform clustering analysis:
```sh
python3 cluster.py <cohort/summary.csv> <group_column1> <group_column2> ... 
```

Perform differentially expression analysis:
```sh
python3 DEA.py <cohort/summary.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... 
```

Perform gene set enrichment analysis:
```sh
python3 GSEA.py <cohort/summary.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... <geneset.gmt>
```

Perform survival analysis:
```sh
python3 survival.py <cohort/summary.csv> <time> <event> <group_column1> <group_column2> ...
```