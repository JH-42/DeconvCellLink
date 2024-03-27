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




DCL_net <- function(expression_data, geneList = NULL, tissueType = NULL, mult = FALSE, mult_tissue = NULL,numCores=2) {
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
      y <- clusterProfiler::enricher(gene = geneList, TERM2GENE = ggs, minGSSize = 3, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
    }
    
    am <- bnpathplot(results = y,
                     exp = expression_data, qvalueCutOff = 0.05,cexLine = 0,
                     R = 100, orgDb = org.Mm.eg.db, nCategory = 50,
                     expRow = "SYMBOL", bypassConverting = T,
                     color = "p.adj", returnNet = TRUE)
    dcl_plot<-am$plot+scale_color_gradient(name = "p.adj",low = "#6DC2C5", high = "#aa0051")
    
    
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
      gs <- melt(df, varnames = c("cellName", "geneID"), 
                 value.name = "expression")
      gs <- gs[, c("geneID", "expression")]
      colnames(gs) <- c("cellName", "geneID")
      row.names(gs) <- gs$SYMBOL
      # symbol <- as.character(gs[, 2])
      # eg_gs <- bitr(symbol, fromType = "SYMBOL", toType = "ENTREZID", 
      #               OrgDb = "org.Mm.eg.db")
      # row.names(gs) <- gs$SYMBOL
      # row.names(eg_gs) <- eg_gs$SYMBOL
      # colnames(eg_gs) <- c("geneID", "ENTREZID")
      # ggs <- merge(gs, eg_gs, by = "geneID")
      # ggs$geneID <- NULL
      # colnames(ggs) <- c("cellName", "geneID")
      # ggs <- ggs %>% distinct(cellName, geneID)
    }
    
    if (is.null(geneList)) {
      return(list(marker = gs, cells_proportion = ep$Proportion))
    } else {
      y <- clusterProfiler::enricher(gene = geneList, TERM2GENE = gs, minGSSize = 3, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
    }  
    y@keytype <- "ENTREZID"
    y@organism <- "Mus musculus"
    
    n1 <- bitr(row.names(expression_data), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    n2 <- expression_data[n1$SYMBOL, ]
    rownames(n2) <- n1$ENTREZID
    
    am <- bnpathplot(results = y,
                     exp = expression_data, qvalueCutOff = 0.05,cexLine = 0,
                     R = 100, orgDb = org.Mm.eg.db, nCategory = 50,
                     expRow = "SYMBOL", bypassConverting = T, returnNet = TRUE)
    
    dcl_plot <- am$plot + scale_color_gradient(name = "p.adj", 
                                                low = "#6DC2C5", high = "#aa0051")
    return(list(enricher = y, marker = gs, bnObject = am, 
                dcl_plot = dcl_plot,
                cells_proportion = combined_prop_df,
                gene_scores = combined_Escore_df,
                tissue=mult_tissue))
  }
}

