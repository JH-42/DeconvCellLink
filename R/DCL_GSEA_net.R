#' Calculate Deconvolution Cell Link for Bulk RNAseq data
#'
#' @param expression_data A data frame containing the expression data (CPM/TPM).
#' @param geneList A list of differentially expressed genes (DEgenes), to be used for SSMD calculation. Default is NULL.
#' @param tissueType The tissue type of the expression data. Choose one of the following: 'Inflammatory', 'Central Nervous System', 'Hematopoietic System', 'Blood'.
#' @param mult Whether to use multi-tissue annotations.
#' @param mult_tissue Which tissue annotations to use.
#' @param numCores How many cores to use.
#' 
#' 
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


DCL_GSEA_net <- function(expression_data, geneList = NULL, tissueType = NULL, mult = FALSE, mult_tissue = NULL, numCores = 2) {
  combined_proportion <- list()
  combined_marker_gene <- list()
  combined_Escore <- list()
  combined_potential_modules <- list()
  
  
  if (!mult) {
    # single tissue
    ep <- SSMD(bulk_data = expression_data, tissue = tissueType)
    gs <- ep$marker_gene
    
    if (is.null(geneList)) {
      return(list(marker = gs, cells_proportion = ep$Proportion))
    } else {
      data_all_sort <- geneList %>% arrange(desc(geneList[1]))
      geneList = data_all_sort$logFC
      names(geneList) <- row.names(data_all_sort)
      symbol_2 <- as.character(rownames(geneList))
      y <- clusterProfiler::GSEA(gene = geneList, TERM2GENE = gs, minGSSize = 3, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = 123)
    }
    
    am <- bnpathplot(results = y,
                     exp = dd, qvalueCutOff = 0.05,
                     R = 100, orgDb = org.Mm.eg.db, nCategory = 50,
                     expRow = "SYMBOL", bypassConverting = T,
                     color = "enrichmentScore", returnNet = TRUE)
    
    return(list(enricher = y, marker = gs, bnObject = am, cells_proportion = combined_prop_df, gene_scores = combined_Escore_df，tissue=tissueType))
  } else {
    # multi tissue
    if (is.null(mult_tissue)) {
      stop("mult_tissue cannot be NULL when mult is TRUE")
    }
    

    # multi_core
    cores <- numCores
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    # foreach cycle
    start_time <- Sys.time()
    combined_results <- foreach::foreach(tissue = mult_tissue, .packages = c("SSMD", "bcv")) %dopar% {
      library(SSMD)
      library(bcv)
      ssmd_result <- SSMD(bulk_data = expression_data, tissue = tissue)
      list(
        Proportion = ssmd_result$Proportion,
        marker_gene = ssmd_result$marker_gene,
        Escore = ssmd_result$Escore
      )
    }
    end_time <- Sys.time()
    cat("SSMD function elapsed time:", end_time - start_time, "\n")
    
    # stop
    stopCluster(cl)
    
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
      data_all_sort <- geneList %>% arrange(desc(geneList[1]))
      geneList = data_all_sort$logFC
      names(geneList) <- row.names(data_all_sort)
      symbol_2 <- as.character(rownames(geneList))
      y <- clusterProfiler::GSEA(gene = geneList, TERM2GENE = gs, minGSSize = 3, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = 123)
    }
    
    am <- bnpathplot(results = y,
                     exp = dd, qvalueCutOff = 0.05,
                     R = 100, orgDb = org.Mm.eg.db, nCategory = 50,
                     expRow = "SYMBOL", bypassConverting = T,
                     color = "enrichmentScore", returnNet = TRUE)
    
    return(list(enricher = y, marker = gs, bnObject = am, 
                cells_proportion = combined_prop_df,
                gene_scores = combined_Escore_df,
                tissue=mult_tissue))
  }
}
