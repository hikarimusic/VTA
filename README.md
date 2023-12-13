# VTA
An variance-tolerant aligner for next-generation sequencing

## Usage
Compile the files:
```sh
make
```
Index the reference genome:
```sh
./index reference_genome.fna
```
Align the sequences:
```sh
./align reference_genome.fna input1.fq input2.fq output.sam
```
Tools:
```sh
./tools reference_genome.fna
```

## Algorithm
* **Index**
    * Count sort
    * Suffix array
    * Burrows-Wheeler Transform
    * FM index
* **Align**
    * Seeding
    * Clustering
    * DP (edit distance)
* **Performance**
    * Index: GRCH38 in 4 hours
    * Align: 1000000 * 200 bp in 2000 seconds
    * Memory: <10 gigabytes