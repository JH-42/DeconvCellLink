#' Title
#'
#' @param expression_data the expression_data input: CPM/TPM
#' @param geneList DEgene list or other gene list， defult=NULL
#' @param tissueType the one of 'Inflammatory', 'Central Nervous System', 'Hematopoietic System','Blood'
#'
#' @return A list containing the following objects:
#'   \item{enrichmentResult}{Enrichment analysis of cell types specific marker.}
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
    print("Error! With out DEgene")
  } else {
    symbol_2 <- as.character(geneList)
    eg_gene <- bitr(symbol_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    geneList <- bitr(geneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    y <- clusterProfiler::enricher(gene = geneList$ENTREZID, TERM2GENE = ggs, minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none")
  }

  y@keytype <- "ENTREZID"
  y@organism <- "Mus musculus"

  n1 <- bitr(row.names(expression_data), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  n2 <- expression_data[n1$SYMBOL, ]
  rownames(n2) <- n1$ENTREZID

  am <- bnpathplot(results = y, exp = n2, expRow = "ENTREZID", orgDb = org.Mm.eg.db, qvalueCutOff = 0.05, adjpCutOff = 0.05, R = 100, seed = 123) # ,interactive=T, strengthPlot = T


  # y <- clusterProfiler::enricher(eg_gene$ENTREZID, TERM2GENE=ggs, minGSSize=1,pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "none")

  return(list(enricher = y, marker = gs, bnObject = am))
}
