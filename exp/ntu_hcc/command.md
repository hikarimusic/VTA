```sh
(.vta) hikari@hikari:~/VTA$ python3 summarize.py hcc_final.csv hcc_dir/ DDX11L2 tpm
[Process TSV] Complete
[Save Result] Complete
(.vta) hikari@hikari:~/VTA$ python3 cluster.py hcc_final/summary.csv Cluster Etiology CTNNB1
[Read Data] Complete
[PCA] Complete
[Filter Genes] 17892
[Hierarchy Cluster] Complete
[Create Heatmap] Complete
[Save Heatmap] Complete
(.vta) hikari@hikari:~/VTA$ python3 DEG.py hcc_final/summary.csv Cluster Cluster1 -- Cluster2 Cluster3
[Read Data] Complete
[Filter Genes] 17892
[Find DEG] Down: 2894 / Up: 5164
[Save Results] Complete
[Chi-square] Group / Etiology: None
[Chi-square] Group / CTNNB1: None
[Create Plots] Complete
(.vta) hikari@hikari:~/VTA$ python3 GSEA.py hcc_final/summary.csv Cluster Cluster1 -- Cluster2 Cluster3 HALLMARK.gmt
[Read Data] Complete
[Filter Genes] 17892
[Load GeneSets] 50
[KS Test] Complete
[Save Results] Complete
[Create Plots] Up: 3 / Down: 10
(.vta) hikari@hikari:~/VTA$ python3 DEG.py hcc_final/summary.csv Cluster Cluster2 -- Cluster3
[Read Data] Complete
[Filter Genes] 17892
[Find DEG] Down: 4154 / Up: 5307
[Save Results] Complete
[Chi-square] Group / Etiology: 2.27e-07
[Chi-square] Group / CTNNB1: 0.366
[Create Plots] Complete
(.vta) hikari@hikari:~/VTA$ python3 GSEA.py hcc_final/summary.csv Cluster Cluster2 -- Cluster3 HALLMARK.gmt
[Read Data] Complete
[Filter Genes] 17892
[Load GeneSets] 50
[KS Test] Complete
[Save Results] Complete
[Create Plots] Up: 2 / Down: 2
```

```sh
pdflatex main.tex
biber main
```