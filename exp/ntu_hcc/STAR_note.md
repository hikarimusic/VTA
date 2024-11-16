```sh
sudo apt-get update
sudo apt-get install fastqc
sudo apt-get install rsem
```
```sh
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
make STAR
```
```sh
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v36.annotation.gtf.gz
```