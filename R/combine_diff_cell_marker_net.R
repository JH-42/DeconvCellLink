#' Build a cell-cell network from a custom cell-type marker table
#'
#' @param expression_data Expression matrix, log(CPM/TPM), with gene symbols as row names.
#' @param geneList A character vector of gene symbols (e.g. differentially expressed genes). Default is NULL.
#' @param cell_markers A marker table with columns 'cellName' and 'geneID'.
#'
#' @return A list containing the following objects:
#'   \item{enricher}{Over-representation result on the cell-type markers.}
#'   \item{marker}{The cell-type marker table used.}
#'   \item{bnObject}{Deconvolution cell link from Bayesian Network.}
#' @export
#'
#'
#'

combine_diff_cell_marker_net <- function(expression_data, geneList = NULL, cell_markers) {
  {
    gs <- cell_markers[, c("cellName", "geneID")]
    colnames(gs) <- c("cellName", "geneID")
    row.names(gs) <- gs$SYMBOL
    symbol <- as.character(gs[, 2])
    eg_gs <- clusterProfiler::bitr(symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    row.names(gs) <- gs$SYMBOL
    row.names(eg_gs) <- eg_gs$SYMBOL
    colnames(eg_gs) <- c("geneID", "ENTREZID")

    ggs <- merge(gs, eg_gs, by = "geneID")
    ggs$geneID <- NULL
    colnames(ggs) <- c("cellName", "geneID")
    ggs <- ggs %>% dplyr::distinct(cellName, geneID)
  }
  if (is.null(geneList)) {
    print("Error! With out DEgene")
  } else {
    symbol_2 <- as.character(geneList)
    eg_gene <- clusterProfiler::bitr(symbol_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    geneList <- clusterProfiler::bitr(geneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    y <- clusterProfiler::enricher(gene = geneList$ENTREZID, TERM2GENE = ggs, minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none")
  }

  y@keytype <- "ENTREZID"
  y@organism <- "Mus musculus"

  n1 <- clusterProfiler::bitr(row.names(expression_data), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  n2 <- expression_data[n1$SYMBOL, ]
  rownames(n2) <- n1$ENTREZID

  am <- CBNplot::bnpathplot(results = y, exp = n2, expRow = "ENTREZID", orgDb = org.Mm.eg.db, qvalueCutOff = 0.05, adjpCutOff = 0.05, R = 100, seed = 123) # ,interactive=T, strengthPlot = T


  # y <- clusterProfiler::enricher(eg_gene$ENTREZID, TERM2GENE=ggs, minGSSize=1,pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "none")

  return(list(enricher = y, marker = gs, bnObject = am))
}
