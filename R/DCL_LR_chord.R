#' Circular chord plot of the deconvolution cell-cell link network
#'
#' Lineages on a circle, arcs are the Bayesian-network cell-cell edges (width = edge
#' count), optionally labelled with their top ligand-receptor pair and filled by a
#' per-lineage score.
#'
#' @param DCL_Object Deconvolution Cell Link object from DCL_net or DCL_GSEA_net.
#' @param node_score Named numeric vector, one value per lineage, to fill nodes. Default NULL.
#' @param lr_table Optional data frame with columns from, to, ligand, receptor (plus a score/cor column for ranking) to label edges. Typically DCL_LR_plot()$results. Default NULL.
#' @param strength_threshold Keep edges with bnlearn strength above this. Default 0.6.
#' @param top_lr Number of ligand-receptor pairs to label per edge. Default 1.
#' @param lineage Collapse subtypes to lineage by dropping a trailing _number. Default TRUE.
#' @param palette Length-3 diverging colours (low, mid, high) for node fill.
#' @param score_limits Length-2 fill limits. Default NULL (symmetric around 0).
#' @param node_size Node point size. Default 8.
#' @param edge_color Edge colour. Default "grey60".
#' @param label_size Text size (mm) for labels. Default 2.1.
#' @param title Plot title. Default NULL.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr %>%
#' @export

DCL_LR_chord <- function(DCL_Object, node_score = NULL, lr_table = NULL,
                         strength_threshold = 0.6, top_lr = 1, lineage = TRUE,
                         palette = c("#2166AC", "white", "#B2182B"), score_limits = NULL,
                         node_size = 8, edge_color = "grey60", label_size = 2.1, title = NULL) {

  lin <- function(x) if (lineage) sub("_[0-9]+$", "", x) else x

  ## cell-cell edges from the BN, filter by strength if present
  bn <- DCL_Object$bnObject
  arcs <- as.data.frame(bn$av$arcs); names(arcs) <- c("from", "to")
  if (!is.null(bn$str) && "strength" %in% names(bn$str)) {
    arcs <- dplyr::inner_join(bn$str, arcs, by = c("from", "to"))
    arcs <- arcs[arcs$strength > strength_threshold, , drop = FALSE]
  }
  edf <- arcs %>% dplyr::mutate(fl = lin(from), tl = lin(to)) %>%
    dplyr::filter(fl != tl) %>% dplyr::count(fl, tl, name = "n")
  if (nrow(edf) == 0) stop("no cell-cell edges left after filtering")

  ## top-N ligand-receptor label per edge
  edf$lab <- ""
  if (!is.null(lr_table) && nrow(lr_table) > 0) {
    lr <- lr_table; lr$fl <- lin(lr$from); lr$tl <- lin(lr$to)
    rank_col <- intersect(c("score", "strength", "comm_score", "cor"), names(lr))[1]
    if (!is.na(rank_col)) lr <- lr[order(-lr[[rank_col]]), ]
    lr <- lr %>% dplyr::filter(fl != tl) %>% dplyr::group_by(fl, tl) %>%
      dplyr::summarise(lab = paste(utils::head(unique(paste0(ligand, "-", receptor)), top_lr),
                                   collapse = "\n"), .groups = "drop")
    edf <- dplyr::left_join(edf, lr, by = c("fl", "tl"), suffix = c("", ".lr"))
    edf$lab <- ifelse(is.na(edf$lab.lr), "", edf$lab.lr); edf$lab.lr <- NULL
  }

  ## lineages evenly on a circle, starting from the top
  lins <- sort(unique(c(edf$fl, edf$tl)))
  ang <- seq(pi/2, pi/2 - 2*pi, length.out = length(lins) + 1)[seq_along(lins)]
  nd <- data.frame(name = lins, x = cos(ang), y = sin(ang), stringsAsFactors = FALSE)
  nd$score <- if (is.null(node_score)) 0 else node_score[nd$name]
  nd$score[is.na(nd$score)] <- 0
  if (is.null(score_limits)) { m <- max(abs(nd$score), 1e-6); score_limits <- c(-m, m) }

  g <- igraph::graph_from_data_frame(edf %>% dplyr::transmute(from = fl, to = tl, n, lab),
                                     vertices = nd, directed = TRUE)
  lay <- ggraph::create_layout(g, layout = "manual", x = nd$x, y = nd$y)

  ggraph::ggraph(lay) +
    ggraph::geom_edge_arc(ggplot2::aes(width = n, label = lab), colour = edge_color,
                          strength = 0.15, alpha = 0.8,
                          arrow = ggplot2::arrow(length = ggplot2::unit(1.5, "mm"), type = "closed"),
                          end_cap = ggraph::circle(5.5, "mm"), start_cap = ggraph::circle(5.5, "mm"),
                          angle_calc = "along", label_size = label_size,
                          label_dodge = ggplot2::unit(1.8, "mm")) +
    ggraph::scale_edge_width(range = c(0.3, 1.5), guide = "none") +
    ggraph::geom_node_point(ggplot2::aes(fill = score), shape = 21, size = node_size,
                            colour = "grey25", stroke = 0.3) +
    ggplot2::scale_fill_gradient2(low = palette[1], mid = palette[2], high = palette[3],
                                  midpoint = 0, limits = score_limits, name = "node score") +
    ggraph::geom_node_text(ggplot2::aes(label = name), size = label_size + 0.4,
                           fontface = "bold", repel = TRUE, max.overlaps = Inf) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(title = title) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(legend.position = "bottom")
}
