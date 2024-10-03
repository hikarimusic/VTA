# VTA
A Complete Library for NGS Data Analysis

## Usage
Compile the files:
```sh
make
```

Install Python packages:
```sh
pip install -r requirements.txt
```

Index the reference genome:
```sh
./index GRCh38.fna
```

Align the sequences:
```sh
./align GRCh38.fna input1.fq input2.fq output.sam
```

Profile the mRNA expression:
```sh
./profile gencode.gtf output.sam output.tsv
```

Analyze the profiles:
```sh
python3 profile.py profile_dir profile_dir/cohort.csv
```

<!-- Tools:
```sh
./tools GRCh38.fna
``` -->

------------------------

*(This library mainly focuses on the analysis of aligned data. If you found the index and align parts too slow, you can use other tools such as [bwa](https://github.com/lh3/bwa) for alignment, and then paste the .sam files here to perform further analysis)*