#' Identify the differentially expressed genes and the corresponding differential expression info p_val, ave_log2FC, percentage of cells expressing that gene and p_val_adj
#'
#' @name get_DEG
#' @param seuratobj the Seurat object of interest
#' @param idents_1 vectors of idents
#' @param idents_2 vectors of idents
#' @param only_pos TRUE or FALSE, indicates whether return positive markers only in each condition in idents_1 compared to each condition in idents_2.
#' @param min_pct a number between 0 and 1. Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations under comparison.
#' @param logfc_threshold a positive number
#' @param p_val_adj_threshold a number between 0 and 1
#' @param test_use Denotes which test to use.  Available options are "wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST" and "DESeq2".
#' @return
#'  DEGinfo: a list of two entries:
#'  DEGinfo$DEG: a list where each entry is a matrix with genes in rows and colomns "p_val", "ave_log2FC", "pct.1", "pct.2" and "p_val_adj"
#'  DEGinfo$DEgenes: a vector of differentially expressed genes' symbols
#' @import nichenetr tidyverse Seurat
#' @export


get_DEG <- function(seuratobj, idents_1, idents_2, only_pos, min_pct, logfc_threshold, p_val_adj_threshold, test_use) {
  Idents(object = seuratobj) = 'Condition'
  DEG =  list()
  for (j in 1:length(idents_2)) {
    for (i in 1:length(idents_1)) {
      seuratmarkers = FindMarkers(seuratobj, ident.1 = idents_1[i], ident.2 = idents_2[j], 
                                  only.pos = only_pos, min.pct = min_pct, logfc.threshold = logfc_threshold, test.use = test_use)
      DEG[[i+(j-1)*length(idents_1)]] = seuratmarkers[which(seuratmarkers[,'p_val_adj'] < p_val_adj_threshold),]
      names(DEG)[i+(j-1)*length(idents_1)] = sprintf("%s_vs_%s", idents_1[i], idents_2[j])
    }
  }
  
  DEgenes = vector()
  for (i in 1:length(DEG)) {
    DEgenes = union(DEgenes, rownames(DEG[[i]]))
  }
  
  DEGinfo = list(DEG, DEgenes)
  names(DEGinfo) = c('DEG', 'DEgenes')
  
  return(DEGinfo)
}








