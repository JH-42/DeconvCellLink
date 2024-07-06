#' Calculate LR in Deconvolution Cell Link
#'
#' @param DCL_Object Deconvolution Cell Link object from DCL_GSEA or DCL_net.
#' @param deg A character of differentially expressed genes (DEgenes), Default is NULL.
#' @param expression_data data frame of expression data
#' @param cor_method correlation method of ligand and receptor, one of "pearson", "kendall" and "spearman"
#' @param cor_threshold the threshold of correlation 
#' 
#' @return A list containing the following objects:
#'   \item{results}{L-R results.}
#'   \item{LR_plot}{L-R plot.}
#'

DCL_LR_plot <- function(DCL_Object=DCL_Object, deg=NULL, expression_data=NULL, cor_method="pearson", cor_threshold=0) {
  
  gene_scores <- data.frame(DCL_Object$gene_scores)
  gene_scores$geneID <- row.names(gene_scores)
  names(gene_scores) <- c("gene_score","geneID")
  
  arcs <- DCL_Object$bnObject$av$arcs
  str_data <- DCL_Object$bnObject$str
  arcs_df <- data.frame(from = arcs[,1], to = arcs[,2])
  cell_interactions <- str_data %>%
    inner_join(arcs_df, by = c("from", "to")) %>%
    filter(strength > 0.6) %>%
    dplyr::select(from, to, strength) %>%
    distinct()
  

  
  cell_genes <- inner_join(DCL_Object$marker, gene_scores, by = "geneID")[,1:2]
  gene_interactions <- LRI_mouse$LRI_curated
  # preprocess LRI
  cell_pairs <- gene_interactions %>%
    separate(LRI, into = c("ligand", "receptor"), sep = ":") %>%
    separate_rows(ligand, sep = "_") %>%
    separate_rows(receptor, sep = "_") %>%
    distinct()

  
  if (length(DCL_Object$tissue) > 1) {
    ssmd_markers <- lapply(sapply(DCL_Object$tissue, switch, 
                                  "Inflammatory" = SSMD::Mouse_Cancer_core_marker,
                                  "Central Nervous System" = SSMD::Mouse_Brain_core_marker,
                                  "Hematopoietic System" = SSMD::Mouse_hematopoietic_core_marker,
                                  "Blood" = SSMD::Mouse_Blood_core_marker), 
                           function(x) x[!is.na(x)])
    ssmd_markers <- unlist(ssmd_markers, recursive = FALSE)
    names(ssmd_markers) <- sub("^.*\\.", "", names(ssmd_markers))
  } else {
    ssmd_markers <- switch(DCL_Object$tissue,
                           "Inflammatory" = SSMD::Mouse_Cancer_core_marker,
                           "Central Nervous System" = SSMD::Mouse_Brain_core_marker,
                           "Hematopoietic System" = SSMD::Mouse_hematopoietic_core_marker,
                           "Blood" = SSMD::Mouse_Blood_core_marker)
  }
  
  results <- data.frame(from = character(),
                        to = character(),
                        strength = numeric(),
                        ligand = character(),
                        receptor = character(),
                        stringsAsFactors = FALSE)
  
  # do all cell pair
  for (i in 1:nrow(cell_pairs)) {
    from_cell <- cell_pairs$from[i]
    to_cell <- cell_pairs$to[i]
    strength <- cell_pairs$strength[i]
    # all gene containing the SSMD geneID
    from_genes <- unique(c(cell_genes$geneID[cell_genes$cellName == from_cell],
                           ssmd_markers[[str_extract(from_cell, "^[[:alpha:]]+((?=_)|$)")]]))
    
    to_genes <- unique(c(cell_genes$geneID[cell_genes$cellName == to_cell],
                         ssmd_markers[[str_extract(to_cell, "^[[:alpha:]]+((?=_)|$)")]]))
    
    
    interactions <- gene_interactions %>%
      filter(ligand %in% from_genes & receptor %in% to_genes)
    
    if (nrow(interactions) > 0) {
      new_rows <- data.frame(from = from_cell,
                             to = to_cell,
                             strength = strength,
                             ligand = interactions$ligand,
                             receptor = interactions$receptor,
                             database = interactions$DATABASE,
                             source = interactions$SOURCE,
                             stringsAsFactors = FALSE)
      results <- bind_rows(results, new_rows)
    }
  }
  
  results <- results %>% distinct(from, to, ligand, receptor, .keep_all = TRUE)
  
  if (!is.null(deg)) {
    deg_genes <- deg
    
    # screening the LR by using deg
    results_in_deg <- results %>%
      filter(ligand %in% deg_genes | receptor %in% deg_genes)
    
    results <- results_in_deg
  }
  
  # Correlation analysis using expression data
  if (!is.null(expression_data)) {
    results <- results[apply(results[, c("ligand", "receptor")], 1, function(x) all(x %in% rownames(expression_data))), ]
    expression_data <- as.matrix(expression_data)
    results_with_cor <- results %>%
      rowwise() %>%
      mutate(cor = cor(expression_data[ligand, ], expression_data[receptor, ], method = cor_method)) %>%
      filter(cor >= cor_threshold) %>%
      ungroup()
    
    results <- results_with_cor
  }
  
  plot_data <- results %>%
    mutate(cell_interaction = paste(from, to, sep = " - "),
           LR_interaction = paste(ligand, receptor, sep = " - ")) %>%
    dplyr::select(cell_interaction, LR_interaction, strength,cor) %>%
    distinct()
  
  cell_interactions <- unique(plot_data$cell_interaction)
  lr_interactions <- unique(plot_data$LR_interaction)
  
  df <- plot_data %>%
    mutate(cell_interaction = factor(cell_interaction, levels = cell_interactions),
           LR_interaction = factor(LR_interaction, levels = lr_interactions))
  
  LR_plot <- ggplot(df, aes(x = cell_interaction, y = LR_interaction, color = strength, size = cor)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size(range = c(1, 5)) +
    labs(x = "Ligand-Receptor Interaction", y = "Cell-Cell Interaction", color = "Strength", size = "Correlation") +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "grey", size = 0.5), 
      panel.grid.minor = element_line(color = "grey", size = 0.25),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "right"
    )
  return(list(results = results, LR_plot = LR_plot))
}
