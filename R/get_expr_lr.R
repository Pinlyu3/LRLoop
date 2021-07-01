#' Identify "expressed" ligand-receptor pairs
#' @name get_expr_lr
#' @param lr_network A ligand-receptor network in matrix with two columns "from" and "to".
#' @param thresh_expr_from A matrix of 0s and 1s, indicates whether each gene is "expressed" in each condition in the sender cells of this step
#' @param thresh_expr_to A matrix of 0s and 1s, indicates whether each gene is "expressed" in each condition in the receiver cells of this step
#' @param conditions A vector of conditions of interest
#' @return
#'  lr_expr: A list with two elements:
#'  lr_expr$eachcondition: A list, each element is a ligand-receptor network matrix with two columns "from" and "to" containing "expressed" ligand-receptor pairs in one condition
#'  lr_expr$bind: A ligand-receptor matrix with two columns "from" and "to" containing "expressed" ligand-receptor pairs in at least one condition
#' @import nichenetr tidyverse Seurat
#' @export



get_expr_lr <- function(lr_network, thresh_expr_from, thresh_expr_to, conditions) {
  
  lr_expr_list = list()
  for (i in 1:length(conditions)) {
    lr_expr_list[[i]] = lr_network[thresh_expr_from[lr_network[,'from'], conditions[i]]==1 & thresh_expr_to[lr_network[,'to'], conditions[i]]==1,]
  }
  names(lr_expr_list) = conditions
  
  lr_expr_bind = vector()
  for (i in 1:length(lr_expr_list)) {
    lr_expr_bind = rbind(lr_expr_bind, lr_expr_list[[i]])
  }
  lr_expr_bind = unique(lr_expr_bind)
  
  lr_expr = list(lr_expr_list, lr_expr_bind)
  names(lr_expr) = c('eachcondition', 'bind')
  
  return(lr_expr)
}


