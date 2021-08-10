#' Calculate the gene detection rate in each condition
#' @name get_pct_expr
#' @param seuratobj A Seurat object 
#' @param conditions A vector of conditions of interest
#' @return
#'  A matrix of the detection rate of each gene in each condition with genes in rows and conditions in columns
#' @import nichenetr tidyverse Seurat
#' @export





get_pct_expr <- function(seuratobj, conditions) {
  library(Seurat)
  library(Matrix)
  data = seuratobj@assays$RNA@data
  pct_expr = matrix(0, nrow = nrow(data), ncol = length(conditions))
  Idents(object = seuratobj) = 'Condition'
  for (i in 1:length(conditions)) {
    subdata = data[,WhichCells(seuratobj, idents = conditions[i])]
    print(dim(subdata))
    pct_expr[,i] = rowSums(subdata >0)/ncol(subdata)
  }
  rownames(pct_expr) = rownames(data)
  colnames(pct_expr) = conditions
  
  return(pct_expr)
}




