# DeconvCellLink: Approach for Infer potential deconvolution cell types interactions from mouse bulk RNA-seq data

## Installation
```
install.packages("remotes")
install.packages("BiocManager")
BiocManager::install("JH-42/DeconvCellLink")
```
## Dependency

```
if (!require(bcv)) install.packages("bcv")
if (!require(SSMD)) BiocManager::install("xiaoyulu95/SSMD")
if (!require(clusterProfiler)) BiocManager::install("clusterProfiler")
if (!require(CBNplot)) BiocManager::install("noriakis/CBNplot")
if (!require(pheatmap)) install.packages("pheatmap")
if (!require(reshape2)) install.packages("reshape2")
if (!require(org.Mm.eg.db)) BiocManager::install("org.Mm.eg.db")
```

## Usage

```
#The exp_dat and geneList should use Gene Symbols.
estimate.proportion <- DCL_net(exp_dat,geneList,tissueType = "Inflammatory")
estimate.proportion$bnObject
pheatmap::pheatmap(estimate.proportion$cells_proportion,scale = "row")
```
## Arguments
* `exp_dat`        A data frame containing the expression data (CPM/TPM).
* `geneList`        A gene list, such as a list of differentially expressed genes (DEgenes).
* `tissueType`        The tissue type of the expression data. Choose one of the following: 'Inflammatory', 'Central Nervous System', 'Hematopoietic System', 'Blood'.

## Value

An object of class is also invisibly returned. This is a list containing
the following components:

* `bnObject`        Deconvolution cell link from Bayesian Network.
* `cells_proportion`        Deconvolution Cell proportion from SSMD.
