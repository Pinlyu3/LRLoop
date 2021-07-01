#' Calculate the LRscore of each "expressed" ligand-receptor pair in each condition
#' @name get_LRscores 
#' @param lr_expr A ligand-receptor network matrix with two columns "from" and "to".
#' @param conditions A vector of conditions of interest
#' @param value.use_from A numeric matrix with genes in rows and conditions in columns, selected as the ligand values to be used in the calculation of LRscores
#' @param value.use_to A numeric matrix with genes in rows and conditions in columns, selected as the receptor values to be used in the calculation of LRscores
#' @param scalar A number, plays the role of a scalar in the calculation of LRscores when the LRscore_method set to 'scsigr' or individual_scale'
#' @param LRscore_method The method of calculating the LRscores, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr' 
#' @return
#'  A matrix with columns "ligand", "receptor" and the condition names that keeps record of the LRscore of each ligand-receptor pair in each condition 
#' @import nichenetr tidyverse Seurat
#' @export



get_LRscores <- function(lr_expr, conditions, value.use_from, value.use_to, scalar, LRscore_method) {
  
  LRscorematrix = matrix(0, nrow = nrow(lr_expr), ncol = (2+length(conditions)))
  colnames(LRscorematrix) = c('from', 'to', conditions)
  for (i in 1:nrow(LRscorematrix)) {
    ligand = lr_expr[i,'from']
    receptor = lr_expr[i,'to']
    lrscores = vector()
    for (j in 1:length(conditions)) {
      lrscores[j] = LRscore(value.use_from[ligand,conditions[j]], value.use_to[receptor,conditions[j]], scalar, LRscore_method)
    }
    LRscorematrix[i,] = c(ligand, receptor, lrscores)
  }
  
  return(LRscorematrix)
}



