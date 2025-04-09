# DeconvCellLink: Approach for Infer potential deconvolution cell types communication from mouse bulk RNA-seq data

## Update

**Version 1.2 (2024-03-17)**

_Added new functions of receptor-ligand interactions_


**Version 1.1 (2024-03-15)**

_Added new functions based on GSEA algorithm_

_Added multi-threading to optimize running speed_


## Installation
```
install.packages("remotes")
install.packages("BiocManager")
BiocManager::install("JH-42/DeconvCellLink")
```

## Usage

```
#The exp_dat and geneList should use Gene Symbols.

#single tissue cell types
DCL_obj <- DCL_net(exp_dat,geneList,tissueType = "Inflammatory")

DCL_obj <- DCL_GSEA_net(expression_data = exp,geneList = geneList,tissueType = NULL,mult = T,
             mult_tissue = "Inflammatory",numCores = 12)

#multiple tissue cell types
DCL_obj <- DCL_GSEA_net(expression_data = exp,geneList = geneList,tissueType = NULL,mult = T,
             mult_tissue = c("Inflammatory","Central Nervous System"),numCores = 12)
DCL_obj$bnObject
```
![DCL Plot](https://github.com/JH-42/DeconvCellLink/blob/main/img/cell-cell.png)

```
pheatmap::pheatmap(DCL_obj$cells_proportion,scale = "row")


#LR_plot
LR_plot<-DCL_GSEA_net_interaction(DCL_obj,deg = deg, expression_data = exp,cor_method = "pearson",cor_threshold = 0)#deg = deg list
LR_plot$LR_plot
```
![LR Plot](https://github.com/JH-42/DeconvCellLink/blob/main/img/LR.png)


## Arguments
* `exp_dat`        A data matrix containing the expression data log(CPM/TPM).
* `geneList`        A gene list, such as a list of differentially expressed genes (DEgenes).
* `tissueType`        The tissue type of the expression data. Choose one of the following: 'Inflammatory', 'Central Nervous System', 'Hematopoietic System', 'Blood'.
*  `mult` Whether to use multi-tissue annotations.
*  `mult_tissue` Which tissue annotations to use.
*  `numCores` How many cores to use.

## Value

An object of class is also invisibly returned. This is a list containing
the following components:

* `bnObject`        Deconvolution cell link from Bayesian Network.
* `cells_proportion`        Deconvolution Cell proportion from SSMD.

## Citation

If you use DeconvCellLink in your research, please cite it as follows:  
**Wang, J.; Zhong, Z.; Luo, H.; Han, Q.; Wu, K.; Jiang, A.; Chen, L.; Gao, Y.; Jiang, Y. Modulation of brain immune microenvironment and cellular dynamics in systemic inflammation. Theranostics 2025, 15 (11), 5153-5171. DOI: 10.7150/thno.107061.**

If you use LR_plot function in your research, please cite this paper as follows:  
**Lagger C, Ursu E, Equey A, Avelar RA, Pisco AO, Tacutu R, de Magalhães JP. scDiffCom: a tool for differential analysis of cell-cell interactions provides a mouse atlas of aging changes in intercellular communication. Nat Aging. 2023 Nov;3(11):1446-1461.** 
