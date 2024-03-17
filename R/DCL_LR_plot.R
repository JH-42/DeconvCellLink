#' Calculate LR in Deconvolution Cell Link
#'
#' @param DCL_Object Deconvolution Cell Link object from DCL_GSEA or DCL_net.
#' @param deg A character of differentially expressed genes (DEgenes), Default is NULL.
#' 
#' 
#' @return A list containing the following objects:
#'   \item{results}{L-R results.}
#'   \item{LR_plot}{L-R plot.}
#'

LR_interaction <- function(DCL_Object=DCL_Object, deg=NULL) {
  
  gene_scores <- data.frame(DCL_Object$gene_scores)
  gene_scores$geneID<-row.names(gene_scores)
  names(gene_scores)<-c("gene_score","geneID")
  

  cell_interactions <- DCL_Object$bnObject$str %>% filter(strength > 0.6)#screening strength >0.6
  cell_genes <- inner_join(DCL_Object$marker, gene_scores, by = "geneID")[,1:2]
  gene_interactions <- LRI_mouse$LRI_curated
  # preprocess LRI
  gene_interactions <- gene_interactions %>%
    separate(LRI, into = c("ligand", "receptor"), sep = ":") %>%
    separate_rows(ligand, sep = "_") %>%
    separate_rows(receptor, sep = "_") %>%
    distinct()
  
  # cell_pairs in cell_interactions 
  cell_pairs <- cell_interactions %>%
    dplyr::select(from, to, strength) %>%
    distinct()
  
  ssmd_markers <- lapply(sapply(DCL_Object$tissue, switch, 
                                "Inflammatory" = SSMD::Mouse_Cancer_core_marker,
                                "Central Nervous System" = SSMD::Mouse_Brain_core_marker,
                                "Hematopoietic System" = SSMD::Mouse_hematopoietic_core_marker,
                                "Blood" = SSMD::Mouse_Blood_core_marker), 
                         function(x) x[!is.na(x)])
  
  ssmd_markers <- unlist(ssmd_markers, recursive = FALSE)
  names(ssmd_markers) <- sub("^.*\\.", "", names(ssmd_markers))
  
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
  
  plot_data <- results %>%
    mutate(cell_interaction = paste(from, to, sep = " - "),
           LR_interaction = paste(ligand, receptor, sep = " - ")) %>%
    dplyr::select(cell_interaction, LR_interaction, strength) %>%
    distinct()
  
  cell_interactions <- unique(plot_data$cell_interaction)
  lr_interactions <- unique(plot_data$LR_interaction)
  
  df <- plot_data %>%
    mutate(cell_interaction = factor(cell_interaction, levels = cell_interactions),
           LR_interaction = factor(LR_interaction, levels = lr_interactions))
  
  LR_plot <- ggplot(df, aes(x = cell_interaction, y = LR_interaction, color = strength)) +
    geom_point(size = 3) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(x = "Ligand-Receptor Interaction", y = "Cell-Cell Interaction", color = "Strength") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(list(results = results, LR_plot = LR_plot))
}
