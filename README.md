# DeconvCellLink

Infer potential cell-cell communication between deconvolution-derived cell types from mouse bulk RNA-seq data.

DeconvCellLink deconvolves a bulk RNA-seq matrix into cell-type proportions and marker genes with SSMD, tests which cell types are associated with a gene list of interest (over-representation with `DCL_net`, or GSEA with `DCL_GSEA_net`), and learns a Bayesian network over those cell types to propose cell-cell links. `DCL_LR_plot` then maps ligand-receptor pairs onto the linked cell types.

## Update

**Version 1.2 (2024-03-17)**: Added ligand-receptor interaction functions.

**Version 1.1 (2024-03-15)**: Added a GSEA-based function, added multi-threading.

## Installation

DeconvCellLink depends on SSMD, which needs **bcv 1.0.2**. bcv has been archived on CRAN,
so `install.packages("bcv")` no longer finds it: install version 1.0.2 from the archive
first. bcv compiles from source, so you also need build tools (**Rtools** on Windows,
Xcode command line tools on macOS).

```r
install.packages(c("remotes", "BiocManager"))

# bcv 1.0.2 from the CRAN archive (required by SSMD)
remotes::install_version("bcv", version = "1.0.2")

# if the line above fails, install the archived source directly:
# install.packages(
#   "https://cran.r-project.org/src/contrib/Archive/bcv/bcv_1.0.2.tar.gz",
#   repos = NULL, type = "source")

BiocManager::install("JH-42/DeconvCellLink")
```

## Input

- **`expression_data`**: a gene-by-sample matrix of log(CPM) or log(TPM) values, with **gene symbols as row names** and samples as columns (mouse).
- **`geneList`**: the gene list of interest. Its form differs between the two functions:
  - `DCL_net` runs over-representation (`clusterProfiler::enricher`), so `geneList` is a **character vector of gene symbols**, e.g. your differentially expressed genes.
  - `DCL_GSEA_net` runs GSEA (`clusterProfiler::GSEA`), so `geneList` is a **data.frame with gene symbols as row names and a ranking metric (e.g. logFC) in the first column**, covering all tested genes. It is sorted internally by that first column.
- **`tissueType`** (single tissue) or **`mult_tissue`** (one or more tissues) selects the SSMD marker reference:

  | value | marker set |
  |---|---|
  | `"Inflammatory"` | inflammation / tumor |
  | `"Central Nervous System"` | brain |
  | `"Hematopoietic System"` | hematopoietic |
  | `"Blood"` | blood |

## Usage

```r
# --- Over-representation from a DEgene set ---
deg <- de_table$gene                                   # character vector of gene symbols
DCL_obj <- DCL_net(expression_data = exp, geneList = deg,
                   tissueType = "Inflammatory")

# --- GSEA from a ranked list ---
ranked <- data.frame(logFC = de_table$logFC,           # rownames = gene, col 1 = ranking metric
                     row.names = de_table$gene)

# one tissue
DCL_obj <- DCL_GSEA_net(expression_data = exp, geneList = ranked,
                        mult = TRUE, mult_tissue = "Inflammatory", numCores = 12)

# several tissues
DCL_obj <- DCL_GSEA_net(expression_data = exp, geneList = ranked,
                        mult = TRUE, mult_tissue = c("Inflammatory", "Central Nervous System"),
                        numCores = 12)

DCL_obj$bnObject          # the cell-cell Bayesian network (edges in $av$arcs)
```
![DCL Plot](https://github.com/JH-42/DeconvCellLink/blob/main/img/cell-cell.png)

```r
# cell-type proportions (cell types x samples)
pheatmap::pheatmap(DCL_obj$cells_proportion, scale = "row")

# ligand-receptor pairs on the linked cell types
LR <- DCL_LR_plot(DCL_obj, deg = deg, expression_data = exp,
                  cor_method = "pearson", cor_threshold = 0)
LR$LR_plot
```
![LR Plot](https://github.com/JH-42/DeconvCellLink/blob/main/img/LR.png)

For a single tissue you can pass either `tissueType = "Inflammatory"` (with `mult = FALSE`, the default) or `mult = TRUE, mult_tissue = "Inflammatory"`, both give the same cell-cell network. Use `mult = TRUE` when you pass more than one tissue.

## Arguments

* `expression_data`: expression matrix, log(CPM/TPM), gene symbols as row names.
* `geneList`: DEgene character vector (`DCL_net`) or ranked data.frame (`DCL_GSEA_net`), see [Input](#input).
* `tissueType`: tissue for the single-tissue run (`mult = FALSE`).
* `mult`: use one or more tissues via `mult_tissue`. Default `FALSE`.
* `mult_tissue`: tissue(s) to use when `mult = TRUE`.
* `numCores`: cores for the parallel SSMD step. Default `2`.
* `seed`: random seed for reproducibility. Default `123`.
* `hub`: optional hub argument passed to the Bayesian network plot.

`DCL_LR_plot(DCL_Object, deg = NULL, expression_data = NULL, cor_method = "pearson", cor_threshold = 0)`

* `DCL_Object`: the object returned by `DCL_net` or `DCL_GSEA_net`.
* `deg`: optional character vector, keeps only ligand-receptor pairs that involve a DEgene.
* `expression_data`: expression matrix, used to correlate each ligand-receptor pair across samples.
* `cor_method`: `"pearson"`, `"kendall"`, or `"spearman"`.
* `cor_threshold`: keep pairs with correlation at or above this value.

## Value

`DCL_net` and `DCL_GSEA_net` return a list with:

* `bnObject`: the cell-cell Bayesian network, edges are in `bnObject$av$arcs`, edge strengths in `bnObject$str`.
* `cells_proportion`: cell-type proportion matrix (cell types x samples), from SSMD.
* `enricher`: the enrichment result (`enrichResult` for `DCL_net`, `gseaResult` for `DCL_GSEA_net`).
* `marker`: cell-type marker table (`cellName`, `geneID`).
* `gene_scores`: SSMD gene scores.
* `dcl_plot`: the network as a ggplot.
* `tissue`: the tissue(s) used.

## Citation

If you use DeconvCellLink in your research, please cite it as follows:  
**Wang J, Zhong Z, Luo H, Han Q, Wu K, Jiang A, Chen L, Gao Y, Jiang Y.  2025. Modulation of brain immune microenvironment and cellular dynamics in systemic inflammation. Theranostics. 15(11):5153-5171.**

If you use LR_plot function in your research, please cite this paper as follows:  
**Lagger C, Ursu E, Equey A, Avelar RA, Pisco AO, Tacutu R, de Magalhães JP. scDiffCom: a tool for differential analysis of cell-cell interactions provides a mouse atlas of aging changes in intercellular communication. Nat Aging. 2023 Nov;3(11):1446-1461.**
