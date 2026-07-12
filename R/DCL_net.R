#' Calculate Deconvolution Cell Link for Bulk RNAseq data
#'
#' @param expression_data Expression matrix, log(CPM/TPM), with gene symbols as row names and samples as columns.
#' @param geneList A character vector of gene symbols (e.g. differentially expressed genes) for over-representation. Default is NULL.
#' @param tissueType Tissue for the single-tissue run (mult = FALSE). One of 'Inflammatory', 'Central Nervous System', 'Hematopoietic System', 'Blood'.
#' @param mult Use one or more tissues via mult_tissue. Default FALSE.
#' @param mult_tissue Tissue(s) to use when mult = TRUE.
#' @param numCores Cores for the parallel SSMD step. Default 2.
#' @param hub Optional hub argument passed to the Bayesian network plot. Default NULL.
#' @param seed Random seed for reproducibility. Default 123.
#'
#' @return A list containing the following objects:
#'   \item{bnObject}{Deconvolution cell link from Bayesian Network.}
#'   \item{cells_proportion}{Deconvolution Cell proportion from SSMD.}
#'
#'
#'
#' @export
#'
#' @importFrom SSMD SSMD
#' @importFrom clusterProfiler enricher
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom dplyr %>%
#' @importFrom foreach %dopar%
#'
#' @references
#' - Jun-Hao Wang#, Zhao-Qian Zhong#, Hai-Hua Luo#, Qi-Zheng Han, Kan Wu, Li Chen, Yanxia Gao*, Yong Jiang*. Comparative Analysis of Brain Inflammation Reveals New Insights into Sepsis Associated Encephalopathy Progression
#'
#'
#'




DCL_net <- function(expression_data, geneList = NULL, tissueType = NULL, mult = FALSE, mult_tissue = NULL,numCores=2,hub=NULL,seed=123) {
  combined_proportion <- list()
  combined_marker_gene <- list()
  combined_Escore <- list()
  combined_potential_modules <- list()
  set.seed(seed)
  expression_data<-as.matrix(expression_data)
  if (!mult) {
    # single tissue
    ep <- SSMD(bulk_data = expression_data, tissue = tissueType)
    # reshape the per-celltype marker_gene list into a 2-col TERM2GENE df for enricher/GSEA
    gs <- data.frame(cellName = rep(names(ep$marker_gene), lengths(ep$marker_gene)),
                     geneID = unlist(ep$marker_gene, use.names = FALSE))

    if (is.null(geneList)) {
      return(list(marker = gs, cells_proportion = ep$Proportion))
    } else {
      y <- clusterProfiler::enricher(gene = geneList, TERM2GENE = gs, minGSSize = 3, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
    }

    am <- CBNplot::bnpathplot(results = y,
                     exp = expression_data, qvalueCutOff = 0.05,cexLine = 0,
                     R = 100, orgDb = org.Mm.eg.db, nCategory = 50,
                     expRow = "SYMBOL", bypassConverting = T,
                     color = "p.adj", returnNet = TRUE)
    dcl_plot <- am$plot + ggplot2::scale_color_viridis_c(option = "C",name = "p.adj")
    
    
    # y <- clusterProfiler::enricher(eg_gene$ENTREZID, TERM2GENE=ggs, minGSSize=1,pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "none")
    
    return(list(enricher = y, marker = gs, bnObject = am, cells_proportion = ep$Proportion, gene_scores = ep$Escore))
  } else {
    # multi tissue
    if (is.null(mult_tissue)) {
      stop("mult_tissue cannot be NULL when mult is TRUE")
    }
    cores <- numCores
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    start_time <- Sys.time()
    doRNG::registerDoRNG(seed)
    combined_results <- foreach::foreach(tissue = mult_tissue,
                                         .packages = c("SSMD", "bcv")) %dopar% {
                                           library(SSMD)
                                           library(bcv)
                                           ssmd_result <- SSMD(bulk_data = expression_data, 
                                                               tissue = tissue)
                                           result <- list(
                                             Proportion = ssmd_result$Proportion,
                                             marker_gene = ssmd_result$marker_gene,
                                             Escore = ssmd_result$Escore)
                                             return(result)
                                         }
    end_time <- Sys.time()
    cat("SSMD function elapsed time:", end_time - start_time,"\n")
    parallel::stopCluster(cl)
    combined_proportion <- lapply(combined_results, function(x) x$Proportion)
    combined_marker_gene <- lapply(combined_results, function(x) x$marker_gene)
    combined_Escore <- lapply(combined_results, function(x) x$Escore)
    combined_prop_df <- do.call(rbind, combined_proportion)
    combined_Escore_df <- do.call(rbind, combined_Escore)
    ep <- list()
    for (sublist in combined_marker_gene) {
      ep <- append(ep, sublist)
    }
    
    df <- as.matrix(do.call(cbind, ep))
    
    {
      gs <- reshape2::melt(df, varnames = c("cellName", "geneID"), 
                 value.name = "expression")
      gs <- gs[, c("geneID", "expression")]
      colnames(gs) <- c("cellName", "geneID")
      row.names(gs) <- gs$SYMBOL
    }
    
    if (is.null(geneList)) {
      return(list(marker = gs, cells_proportion = ep$Proportion))
    } else {
      y <- clusterProfiler::enricher(gene = geneList, TERM2GENE = gs, minGSSize = 3, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
    }  
    y@keytype <- "ENTREZID"
    y@organism <- "Mus musculus"
    
    n1 <- clusterProfiler::bitr(row.names(expression_data), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    n2 <- expression_data[n1$SYMBOL, ]
    rownames(n2) <- n1$ENTREZID

    am <- CBNplot::bnpathplot(results = y,
                     exp = expression_data, qvalueCutOff = 0.05,cexLine = 0,
                     R = 100, orgDb = org.Mm.eg.db, nCategory = 50,
                     expRow = "SYMBOL", bypassConverting = T, returnNet = TRUE,hub=hub)

    dcl_plot <- am$plot + ggplot2::scale_color_manual(colorRampPalette(c("#F5FB67", "#9E2EA4", "#0D0887"))(100))
    return(list(enricher = y, marker = gs, bnObject = am, 
                dcl_plot = dcl_plot,
                cells_proportion = combined_prop_df,
                gene_scores = combined_Escore_df,
                tissue=mult_tissue))
  }
}

