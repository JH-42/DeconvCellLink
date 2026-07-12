#' Compare functional enrichment across deconvolution cell types
#'
#' @param DCL_Object The object returned by DCL_net or DCL_GSEA_net.
#' @param celltypes Cell types to compare, matched by name. E.g. "Monocyte".
#' @param enrichType One of 'GO', 'Reactome', 'KEGG'.
#' @param showCategory Number of categories to show per cell type. Default 5.
#' @param pvalueCutoff P-value cutoff for enrichment. Default 1.
#' @param pAdjustMethod Multiple-testing adjustment method. Default "BH".
#' @param minGSSize Minimum gene set size. Default 1.
#' @param qvalueCutoff Q-value cutoff for enrichment. Default 1.
#'
#' @return A list containing the following objects:
#'   \item{enrichmentResult}{compareCluster result across the chosen cell types.}
#'   \item{comparePlot}{Dot plot of the compared enrichment.}
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
    dplyr::group_by(cellName) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    tidyr::pivot_wider(names_from = cellName, values_from = geneID) %>%
    dplyr::select(-row)

  qq <- as.data.frame(df_pivoted)

  col_names <- colnames(qq)

  q_col <- qq[grep(celltypes, col_names)]

  q_list <- as.list(q_col)

  clean_list <- function(lst) {
    lapply(lst, function(x) {
      x <- unique(na.omit(x))
      clusterProfiler::bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    })
  }

  qqlist <- clean_list(q_list)

  qalist <- lapply(qqlist, function(x) as.character(x[["ENTREZID"]]))

  if (enrichType == "GO") {
    qcomp <- clusterProfiler::compareCluster(qalist, fun = clusterProfiler::enrichGO, OrgDb = org.Mm.eg.db, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = 1, qvalueCutoff = qvalueCutoff)
  } else if (enrichType == "Reactome") {
    qcomp <- clusterProfiler::compareCluster(qalist, fun = ReactomePA::enrichPathway, organism = "mouse", pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = 1, qvalueCutoff = qvalueCutoff)
  } else if (enrichType == "KEGG") {
    qcomp <- clusterProfiler::compareCluster(qalist, fun = clusterProfiler::enrichKEGG, organism = "mmu", pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = 1, qvalueCutoff = qvalueCutoff)
  }


  cell_comp_plot <- enrichplot::dotplot(qcomp, showCategory = showCategory, includeAll = TRUE) + ggplot2::scale_color_continuous(low = "#9400D3", high = "#DB7093") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  return(list(enrichmentResult = qcomp, comparePlot = cell_comp_plot))
}
