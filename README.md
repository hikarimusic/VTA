# VTA
A Complete Library for NGS Data Analysis

## Usage
Compile the files:
```sh
make
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

<!-- Tools:
```sh
./tools GRCh38.fna
``` -->

## Algorithm
* **Index**
    * Count sort
    * Suffix array
    * Burrows-Wheeler Transform
    * FM index
* **Align**
    * Seeding
    * Clustering
    * DP
* **Performance**
    * Index: GRCH38 in 4 hours
    * Align: 1000000 * 200 bp in 2000 seconds
    * Memory: <10 gigabytes