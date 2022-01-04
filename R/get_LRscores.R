#' Calculate the LRscore of each "expressed" ligand-receptor pair in each condition.  Only run after running "PrepareBasics".
#' @name get_LRscores 
#' @param lr_expr A ligand-receptor network matrix with two columns "from" and "to".
#' @param conditions A vector of conditions of interest
#' @param value.use_from A numeric matrix with genes in rows and conditions in columns, selected as the ligand values to be used in the calculation of LRscores
#' @param value.use_to A numeric matrix with genes in rows and conditions in columns, selected as the receptor values to be used in the calculation of LRscores
#' @param scalar A number, plays the role of a scalar in the calculation of LRscores when the LRscore_method set to 'scsigr' or individual_scale'
#' @param LRscore_method The method of calculating the LRscores, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr' 
#' @param thresh_expr_from A numeric matrix of 0s and 1s with genes in rows and conditions in columns (with the same column names as value.use_from) that define whether each gene is expressed in each condition (0-NO, 1-YES) in the sender cells.
#' @param thresh_expr_to A numeric matrix of 0s and 1s with genes in rows and conditions in columns (with the same column names as value.use_to) that define whether each gene is expressed in each condition (0-NO, 1-YES) in the receiver cells.
#' @param LRL_eachcondition List of LRloop networks in each condition.
#' @param LRL_filter "L1R1", "L2R2" or "none".  "L1R1" ("L2R2"): The LR pairs are from ct1 (ct2) to ct2 (ct1) and apply the LRloop filter in each condition, that is, in each condition, an LRscore will be set to 0 if the LR pair does not form any LRloop in that conditon. "none": No LRloop filter applied. 
#' @param thresh_expr_cut TRUE or FALSE.  If TRUE, the LRscore of ligand-receptor pairs for which either the ligand or the recepor is not expressed as defined in thresh_expr_from and thresh_expr_to, respectively, is set to 0.  Default is FALSE.
#' @return
#'  A matrix with columns "ligand", "receptor" and the condition names that keeps record of the LRscore of each ligand-receptor pair in each condition 
#' @import nichenetr tidyverse Seurat
#' @export


get_LRscores <- function(lr_expr, conditions, value.use_from, value.use_to, scalar, LRscore_method, thresh_expr_from, thresh_expr_to,
                         LRL_eachcondition, LRL_filter = 'none', thresh_expr_cut = FALSE) {
  
  LRscorematrix = matrix(0, nrow = nrow(lr_expr), ncol = (2+length(conditions)))
  colnames(LRscorematrix) = c('from', 'to', conditions)
  for (i in 1:nrow(LRscorematrix)) {
    ligand = lr_expr[i,'from']
    receptor = lr_expr[i,'to']
    lrscores = vector()
    for (j in 1:length(conditions)) {
      if (thresh_expr_cut == FALSE) {
        lrscores[j] = LRscore(value.use_from[ligand,conditions[j]], value.use_to[receptor,conditions[j]], scalar, LRscore_method)
      } else if (thresh_expr_cut == TRUE) {
        lrscores[j] = LRscore(value.use_from[ligand,conditions[j]], value.use_to[receptor,conditions[j]], scalar, LRscore_method)*thresh_expr_from[ligand,conditions[j]]*thresh_expr_to[receptor,conditions[j]]
      }
    }
    LRscorematrix[i,] = c(ligand, receptor, lrscores)
  }
  
  if (LRL_filter == "L1R1") {
    for (j in 1:length(conditions)) {
      cond = conditions[j]
      LRset = unique(LRL_eachcondition[[cond]][['R1->L2 & R2->L1']][,c("L1", "R1")])
      for (i in 1:nrow(LRscorematrix)) {
        idx = which(LRset[,'L1']==LRscorematrix[i,'from'] & LRset[,'R1']==LRscorematrix[i,'to'])
        LRscorematrix[i,cond] = as.numeric(LRscorematrix[i,cond])*(length(idx)>0)
      }
    }
  } else if (LRL_filter == "L2R2") {
    for (j in 1:length(conditions)) {
      cond = conditions[j]
      LRset = unique(LRL_eachcondition[[cond]][['R1->L2 & R2->L1']][,c("L2", "R2")])
      for (i in 1:nrow(LRscorematrix)) {
        idx = which(LRset[,'L2']==LRscorematrix[i,'from'] & LRset[,'R2']==LRscorematrix[i,'to'])
        LRscorematrix[i,cond] = as.numeric(LRscorematrix[i,cond])*(length(idx)>0)
      }
    }
  }
  
  
  
  return(LRscorematrix)
}



