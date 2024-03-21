#' Title
#'
#' @param celltypeObject the object from celltype or celltype_combine
#' @param celltypes choose your favorite cell types E.g: "Monocyte"
#' @param enrichType the one of 'GO', 'Reactome', 'KEGG'
#'
#' @return A list containing the following objects:
#'   \item{bnObject}{Deconvolution cell link from Bayesian Network.}
#'   \item{cells_proportion}{Deconvolution Cell proportion from SSMD.}
#' @export
#'
#'
#'


# options(clusterProfiler.download.method = "wininet")#windows
# options(clusterProfiler.download.method = "auto")#windows
# options(clusterProfiler.download.method = "wget")#windows


compareEnrichment <- function(DCL_Object=DCL_Object, celltypes, enrichType, showCategory = 5, pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 1, qvalueCutoff = 1) {
  df <- DCL_Object$marker
  # Pivot data frame
  df_pivoted <- df %>%
    group_by(cellName) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = cellName, values_from = geneID) %>%
    dplyr::select(-row)

  qq <- as.data.frame(df_pivoted)

  col_names <- colnames(qq)

  q_col <- qq[grep(celltypes, col_names)]

  q_list <- as.list(q_col)

  clean_list <- function(lst) {
    lapply(lst, function(x) unique(na.omit(x)))
    lapply(lst, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
  }

  qqlist <- clean_list(q_list)

  qalist <- lapply(qqlist, function(x) as.character(x[["ENTREZID"]]))

  if (enrichType == "GO") {
    qcomp <- compareCluster(qalist, fun = enrichGO, OrgDb = org.Mm.eg.db, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = 1, qvalueCutoff = qvalueCutoff)
  } else if (enrichType == "Reactome") {
    qcomp <- compareCluster(qalist, fun = "enrichPathway", organism = "mouse", pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = 1, qvalueCutoff = qvalueCutoff)
  } else if (enrichType == "KEGG") {
    qcomp <- compareCluster(qalist, fun = "enrichKEGG", organism = "mmu", pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = 1, qvalueCutoff = qvalueCutoff)
  }


  cell_comp_plot <- dotplot(qcomp, showCategory = showCategory, includeAll = TRUE) + scale_color_continuous(low = "#9400D3", high = "#DB7093") +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(list(enrichmentResult = qcomp, comparePlot = cell_comp_plot))
}
