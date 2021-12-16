# shiny_enrichment_analysis
shiny_enrichment_analysis projet 


**Install biomaRT package :** 

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")
```


**Install pathview package :** 

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("pathview")
```

**Install clusterProfiler package :**

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")
```
**Install org.Dr.eg.db package :**
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Dr.eg.db")
```
***P value ajustement ****

```
p.adjust(p, method = p.adjust.methods, n = length(p))
p.adjust.methods
c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
"fdr", "none")
```
à lire :

https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html
https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html
https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
📖 Introduction | Biomedical Knowledge Mining using GOSemSim and clu...
Biomedical knowledge mining using GOSemSim and clusterProfiler.

https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
http://www.jybaudot.fr/Inferentielle/kolmogorov.html

https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html#ora-algorithm 

```
install.packages("golem")
```
