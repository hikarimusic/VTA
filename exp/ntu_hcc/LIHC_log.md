```sh
(.amaterasu) hikari@hikari:~/AMATERASU$ python3 cluster.py LIHC/summary_cluster.csv Cluster Tumor Stage Race
[Read Data] Complete
[PCA] Complete
[Filter Genes] 15896
[Hierarchy Cluster] Complete
[Create Heatmap] Complete
[Save Heatmap] Complete
```

```sh
(.amaterasu) hikari@hikari:~/AMATERASU$ python3 DEA.py LIHC/summary_cluster.csv Cluster Cluster1 -- Cluster2
[Read Data] Complete
[Filter Genes] 15896
[Find DEG] Up: 588 / Down: 1844
[Save Genes] Complete
[Volcano Plot] Complete
[Strip Plot] Complete
[Heatmap] Complete
(.amaterasu) hikari@hikari:~/AMATERASU$ python3 DEA.py LIHC/summary_cluster.csv Cluster Cluster1 -- Cluster3
[Read Data] Complete
[Filter Genes] 15896
[Find DEG] Up: 819 / Down: 2025
[Save Genes] Complete
[Volcano Plot] Complete
[Strip Plot] Complete
[Heatmap] Complete
(.amaterasu) hikari@hikari:~/AMATERASU$ python3 DEA.py LIHC/summary_cluster.csv Cluster Cluster2 -- Cluster3
[Read Data] Complete
[Filter Genes] 15896
[Find DEG] Up: 859 / Down: 737
[Save Genes] Complete
[Volcano Plot] Complete
[Strip Plot] Complete
[Heatmap] Complete
```

```sh
(.amaterasu) hikari@hikari:~/AMATERASU$ python3 GSEA.py LIHC/summary_cluster.csv Cluster Cluster1 -- Cluster2 Cluster3 H
ALLMARK.gmt
[Read Data] Complete
[Filter Genes] 15896
[Load GeneSets] 50
[KS Test] Up: 4 / Down: 3
[Save GeneSets] Complete
[GSEA Plots] Complete
[Bar Plots] Complete
(.amaterasu) hikari@hikari:~/AMATERASU$ python3 GSEA.py LIHC/summary_cluster.csv Cluster Cluster2 -- Cluster1 Cluster3 H
ALLMARK.gmt
[Read Data] Complete
[Filter Genes] 15896
[Load GeneSets] 50
[KS Test] Up: 3 / Down: 7
[Save GeneSets] Complete
[GSEA Plots] Complete
[Bar Plots] Complete
(.amaterasu) hikari@hikari:~/AMATERASU$ python3 GSEA.py LIHC/summary_cluster.csv Cluster Cluster3 -- Cluster1 Cluster2 H
ALLMARK.gmt
[Read Data] Complete
[Filter Genes] 15896
[Load GeneSets] 50
[KS Test] Up: 8 / Down: 5
[Save GeneSets] Complete
[GSEA Plots] Complete
[Bar Plots] Complete
```

```sh
(.amaterasu) hikari@hikari:~/AMATERASU$ python3 survival.py LIHC/summary_cluster.csv Days Event Cluster Tumor Race PCLAF
 -all
[Read Data] Complete
[Analyze Variables] Tumor: 2.710e-08
[Analyze Variables] Stage: 3.479e-06
[Analyze Variables] Race: 3.100e-01
[Analyze Variables] Gender: 2.569e-01
[Analyze Variables] Cluster: 8.045e-02
[Analyze Variables] PCLAF: 8.818e-06
[Analyze Genes] Significant genes: 172
``` 