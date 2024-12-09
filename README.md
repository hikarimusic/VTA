# AMATERASU
[![License](https://img.shields.io/github/license/hikarimusic/AMATERASU)](https://github.com/hikarimusic/AMATERASU/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/633322949.svg)](https://doi.org/10.5281/zenodo.14337457)
[![Stars](https://img.shields.io/github/stars/hikarimusic/AMATERASU)](https://github.com/hikarimusic/AMATERASU)

Automated Mapping And Transcriptome Expression Research Analysis Suite

## Setup

Clone the repository:

```sh
git clone https://github.com/hikarimusic/AMATERASU.git
```

Make a virtual environment:

```sh
python3 -m venv .amaterasu
source .amaterasu/bin/activate
```

Install python packages:

```sh
cd AMATERASU/
pip install -r requirements.txt
```

Download and index genome files with STAR:

```sh
chmod +x setup_STAR.sh process_STAR.sh
./setup_STAR.sh
```

## Raw Data Processing

Start AMATERASU:

```sh
source .amaterasu/bin/activate
cd AMATERASU/
```

Align and profile RNA-seq data with STAR:

```sh
./process_STAR.sh <input_dir/> <output_dir/>
```

The `<input_dir/>` should contain sequencing raw data. Both pair-end and single-end data are acceptable. Example of `<input_dir/>`:

```sh
my_sequences/
├── sample1_1.fastq # paired-end data
├── sample1_2.fastq
├── sample2.fastq # single-end data
└── sample3/
    ├── sample3_1.fastq # Files can be in subfolders
    └── sample3_2.fastq
```

The `<output_dir/>` will contain expression profiles of the samples. Example of `<output_dir/>`:

```sh
my_profiles/
├── sample1.tsv
├── sample2.tsv
└── sample3.tsv
```

<!-- ```sh
make
./index <GRCh38.fna>
./align <GRCh38.fna> <input1.fq> <input2.fq> <output.sam>
./profile <gencode.gtf> <output.sam> <output.tsv>
``` -->

## Cohort Summary

Start AMATERASU:

```sh
source .amaterasu/bin/activate
cd AMATERASU/
```

Summarize the profiles based on your cohort:

```sh
python3 summarize.py <cohort.csv> <profile_dir/> <value_type>
```

The first column of the cohort should be the sample ID. Example of `<cohort.csv>` (or `<cohort.tsv>`):

| sample_id | gender | race | tumor_stage | 
| :- | :- | :- | :- |
| TCGA_DD_AAVP | male | asian | Stage I |
| TCGA_WX_AA46 | male | white | Stage II |
| TCGA_BD_A3EP | female | black | Stage I |

The `<profile_dir/>` should contain profiles of those samples:

```sh
my_profiles/
├── TCGA_DD_AAVP.tsv
├── TCGA_WX_AA46.tsv
└── TCGA_BD_A3EP.tsv
```

`<value_type>` specifies the expression value to be used. Example of a profile `.tsv` file (in this case we can use `tpm_unstranded` as the value type):

```
# gene-model: GENCODE v36
gene_id	gene_name	gene_type	unstranded	stranded_first	stranded_second	tpm_unstranded	fpkm_unstranded	fpkm_uq_unstranded
N_unmapped			8120268	8120268	8120268			
N_multimapping			4402480	4402480	4402480			
N_noFeature			1446120	24573336	24781939			
N_ambiguous			5551010	1476463	1462076			
ENSG00000000003.15	TSPAN6	protein_coding	3231	1634	1597	44.6685	16.8557	23.8956
ENSG00000000005.6	TNMD	protein_coding	0	0	0	0.0000	0.0000	0.0000
ENSG00000000419.13	DPM1	protein_coding	1083	556	528	56.2676	21.2327	30.1006
```

After running the command, a folder `<cohort/>` will be created. The original cohort with gene expression values will be generated as `<cohort/summary.csv>`, which can be used in the following analysis. Example of `<cohort/summary.csv>`:

| sample_id | ... | TSPAN6 | TNMD | ... |
| :- | :- | :- | :- | :- |
| TCGA_DD_AAVP | ... | 136.4 | 0.00 | ... |
| TCGA_WX_AA46 | ... | 45.3 | 0.19 | ... |
| TCGA_BD_A3EP | ... | 60.1 | 0.04 | ... |

## Clustering Analysis

Start AMATERASU:

```sh
source .amaterasu/bin/activate
cd AMATERASU/
```

Perform clustering analysis with columns in interest:

```sh
python3 cluster.py <cohort/summary.csv> <group_column1> <group_column2> ... 
```

Each of the `<group_column?>` will be labeled in the clustering results.

PCA plots corresponding to each column will be generated as `<cohort/cluster_PCA_....png>`. Example:

<img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/cluster_PCA.png" height=300>

Heatmaps with hierarchy clustering will be generated as `<cohort/cluster_heatmap_....png>`. Example:

<img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/cluster_heatmap.png" height=600>

If `[-n]` is specified in the command, samples will be labeled into n clusters and `cohort/summary_cluster.csv` will be generated. You can use it instead of `cohort/summary.csv` in the following analysis to access the `Cluster` column. Example of `cohort/summary_cluster.csv`:

| sample_id | gender | ... | Cluster | 
| :- | :- | :- | :- |
| TCGA_DD_AAVP | male | ... | Cluster1 |
| TCGA_WX_AA46 | male | ... | Cluster2 |
| TCGA_BD_A3EP | female | ... | Cluster3 |

All the configuration such as figure size and format can be set in the head of `cluster.py`.

## Differential Expression Analysis

Start AMATERASU:

```sh
source .amaterasu/bin/activate
cd AMATERASU/
```

Perform differential expression analysis between two groups:

```sh
python3 DEA.py <cohort/summary.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... 
```

The first group will contain samples with `<group_column>` equal to `<group_a?>`, and the second group will contain samples with `<group_column>` equal to `<group_b?>`.

Volcano plots and heatmaps comparing the two groups will be generated as `<cohort/DEA_volcano_....png>` and `<cohort/DEA_heatmap_....png>`. Example:

<img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/DEA_volcano.png" height=300><img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/DEA_heatmap.png" height=300>

Strip plots of the differential expression genes will be generated as `<cohort/DEA_strip_....png>`. Example:

<img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/DEA_strip.png" height=600>

The differential expression genes will be summarized as `cohort/DEA_genes_....csv`. Example:

| gene | log2_fold_change | p_value | adjusted_pvalue |
| :- | :- | :- | :- |
| SAR1B | -1.098911206284902 | 7.924706080956533e-44 | 1.2597112786288506e-39 |
| ACSM2A | -2.25593421738425 | 1.080155659882675e-40 | 8.585077184747501e-37 |
| ACSM2B | -2.023553155838623 | 2.481105941588788e-40 | 1.3146553349165125e-36 |

All the configuration such as figure size and format can be set in the head of `DEA.py`.

## Gene Set Enrichment Analysis

Start AMATERASU:

```sh
source .amaterasu/bin/activate
cd AMATERASU/
```

Perform gene set enrichment analysis between two groups based on the predefined gene set:

```sh
python3 GSEA.py <cohort/summary.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... <geneset.gmt>
```

The first group will contain samples with `<group_column>` equal to `<group_a?>`, and the second group will contain samples with `<group_column>` equal to `<group_b?>`. Predefined gene sets and their genes should be specified in `<geneset.gmt>`. Example:

```
HALLMARK_ADIPOGENESIS    <url>    ABCA1    ABCB8    ACAA2    ...
HALLMARK_ALLOGRAFT_REJECTION    <url>    AARS1    ABCE1    ABI1    ...
HALLMARK_ANDROGEN_RESPONSE    <url>    ABCC4    ABHD2    ACSL3    ...
```

GSEA plots of significant gene sets will be generated in `<cohort/GSEA_gsea_.../>` as `<up_....png>` or `<down_....png>`. Example:

<img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/GSEA_gsea_up.png" height=300><img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/GSEA_gsea_down.png" height=300>

Bar plots showing leading genes of significant gene sets will be generated in `<cohort/GSEA_bar_.../>` as `<up_....png>` or `<down_....png>`. Example:

<img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/GSEA_bar_up.png" height=600><img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/GSEA_bar_down.png" height=600>

The gene sets will be summarized as `cohort/GSEA_genesets_....csv`. Example:

| gene_set | enrichment_score | position | p_value | adjusted_pvalue |
| :- | :- | :- | :- | :- |
| HALLMARK_XENOBIOTIC_METABOLISM | -0.350 | 0.639 | 2.31e-21 | 1.15e-19 |
| HALLMARK_BILE_ACID_METABOLISM | -0.457 | 0.657 | 1.27e-20 | 3.19e-19 |
| HALLMARK_MYC_TARGETS_V1 | 0.336 | 0.240 | 2.92e-20 | 4.86e-19 |

All the configuration such as figure size and format can be set in the head of `GSEA.py`.

## Survival Analysis

Start AMATERASU:

```sh
source .amaterasu/bin/activate
cd AMATERASU/
```

Perform survival analysis to compare among different groups:
```sh
python3 survival.py <cohort/summary.csv> <time> <event> <group_column1> <group_column2> ...
```

`<time>` should record the event time, and `<event>` should record the event status (0 for alive, 1 for dead). `<group_column?>` can be group columns or gene names to be compared.

Kaplan–Meier plots corresponding to each column will be generated as `<cohort/survival_KM_....png>`. Example:

<img src="https://github.com/hikarimusic/AMATERASU/raw/main/assets/survival_KM.png" height=300>

The survival differences will be summarized as `<cohort/survival_variables.csv>`. Example:

| variable | type | threshold | df | p_value |
| :- | :- | :- | :- | :- |
| Tumor | categorical |  | 3 | 2.71e-08 |
| PCLAF | numerical | 2.29 | 1 | 8.82e-06 |

If `[-all]` is specified, all genes will be analyzed and summarized as `<cohort/survival_genes.csv>` (_instead of finding optimal thresholds, the medians will be used_). Example:

| gene | threshold | df | p_value | adjusted_pvalue |
| :- | :- | :- | :- | :- |
| LYPLA2 | 99.9 | 1 | 2.89e-06 | 0.00288 |
| MRTO4 | 30.6 | 1 | 7.68e-06 | 0.00378 |
| YBX1 | 184.2 | 1 | 1.14e-05 | 0.00378 |

All the configuration such as figure size and format can be set in the head of `survival.py`.

## Others

* **Contact**: hikarimusic.tm@gmail.com
* **Citation**: Tung, Yeu-Guang (2024). AMATERASU: Automated Mapping And Transcriptome Expression Research Analysis Suite. Zenodo. https://doi.org/10.5281/zenodo.14337458