#' Calculate Deconvolution Cell Link for Bulk RNAseq data
#'
#' @param expression_data A data frame containing the expression data (CPM/TPM).
#' @param geneList A gene list, such as a list of differentially expressed genes (DEgenes), to be used for SSMD calculation. Default is NULL.
#' @param tissueType The tissue type of the expression data. Choose one of the following: 'Inflammatory', 'Central Nervous System', 'Hematopoietic System', 'Blood'.
#'
#' @return A list containing the following objects:
#'   \item{bnObject}{Deconvolution cell link from Bayesian Network.}
#'   \item{cells_proportion}{Deconvolution Cell proportion from SSMD.}
#'
#'
#'
#' @export
#'
#' @import BiocManager
#' @import remotes
#' @importFrom SSMD SSMD
#' @importFrom clusterProfiler enricher
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @description
#' This package provides functionality for [brief description of what your package does].
#'
#' To install the required packages, please make sure you have BiocManager and remotes installed. If not, you can install them using the following code:
#' \preformatted{
#' if (!require(BiocManager))
#'   install.packages("BiocManager")
#'
#' if (!require(remotes))
#'   install.packages("remotes")
#' }
#'
#' To install this package, you can use the following code:
#' \preformatted{
#' BiocManager::install("JH-42/DevconCellLink")
#' }
#'
#' @references
#' - Jun-Hao Wang#, Zhao-Qian Zhong#, Hai-Hua Luo#, Qi-Zheng Han, Kan Wu, Li Chen, Yanxia Gao, Yong Jiang*. Comparative Analysis of Brain Inflammation Reveals New Insights into Sepsis Associated Encephalopathy Progression
#'
#'
#'

.onLoad <- function(libname, pkgname) {
names(Mouse_Brain_core_marker)<-c("ependymal","micro_glial","oligodendrocyte",
                                  "stromal_like_cell","endothelial","Schwann",
                                  "glial_cell","neuron","astrocyte")

colnames(Mouse_Brain_labeling_matrix)<-c("astrocyte","endothelial","Ependymal",
                                         "stromal_like_cell","oligodendrocyte",
                                         "Microglial","glial_cell","neuron","Schwann")
}

DCL_net <- function(expression_data, geneList = NULL, tissueType, hub = 3) {
  if (tissueType == "Inflammatory") {
    ep <- SSMD(bulk_data = expression_data, tissue = "Inflammatory")
  }
  if (tissueType == "Central Nervous System") {
    ep <- SSMD(bulk_data = expression_data, tissue = "Central Nervous System")
  }
  if (tissueType == "Hematopoietic System") {
    ep <- SSMD(bulk_data = expression_data, tissue = "Hematopoietic System")
  }
  if (tissueType == "Blood") {
    ep <- SSMD(bulk_data = expression_data, tissue = "Blood")
  }


  df <- as.matrix(do.call(cbind, ep$marker_gene))
  
  {
    gs <- melt(df, varnames = c("cellName", "geneID"), value.name = "expression")
    gs <- gs[, c("geneID", "expression")]
    colnames(gs) <- c("cellName", "geneID")
    row.names(gs) <- gs$SYMBOL
    symbol <- as.character(gs[, 2])
    eg_gs <- bitr(symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    row.names(gs) <- gs$SYMBOL
    row.names(eg_gs) <- eg_gs$SYMBOL
    colnames(eg_gs) <- c("geneID", "ENTREZID")
    
    ggs <- merge(gs, eg_gs, by = "geneID")
    ggs$geneID <- NULL
    colnames(ggs) <- c("cellName", "geneID")
    ggs <- ggs %>% distinct(cellName, geneID)
  }
  if (is.null(geneList)) {
    return(list(marker = gs, cells_proportion = ep$Proportion))
  } else {
    symbol_2 <- as.character(geneList)
    eg_gene <- bitr(symbol_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    geneList <- bitr(geneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    y <- clusterProfiler::enricher(gene = geneList$ENTREZID, TERM2GENE = ggs, minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
  }
  
  y@keytype <- "ENTREZID"
  y@organism <- "Mus musculus"
  
  n1 <- bitr(row.names(expression_data), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  n2 <- expression_data[n1$SYMBOL, ]
  rownames(n2) <- n1$ENTREZID
  
  am <- bnpathplot(results = y, exp = n2, expRow = "ENTREZID", orgDb = org.Mm.eg.db, qvalueCutOff = 0.05, adjpCutOff = 0.05, R = 100, seed = 123,hub=hub) # ,interactive=T, strengthPlot = T
  

  # y <- clusterProfiler::enricher(eg_gene$ENTREZID, TERM2GENE=ggs, minGSSize=1,pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "none")

  return(list(enricher = y, marker = gs, bnObject = am, cells_proportion = ep$Proportion))
}
