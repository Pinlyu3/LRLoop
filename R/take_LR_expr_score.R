#' Create a dataframe for the expr_scores of L1R1 or L2R2 pairs across all conditions stored in the list resulted from the function 'LRL_info_collection'
#' @name take_LR_expr_score
#' @param LRloop_info: A list resulted from the function 'LRL_info_collection' 
#' @param LRpair: 'L1R1' or "L2R2'
#' @import nichenetr tidyverse Seurat
#' @export



take_LR_expr_score <- function(LRloop_info, LRpair) {
  if (LRpair == 'L1R1') {
    scorematrix = LRloop_info$`L1R1 expr_score`
    pairnames = sprintf("%s_%s", scorematrix[,'L1'], scorematrix[,'R1'])
    scorevalues = scorematrix[,c(3:(ncol(scorematrix)-2))]
    class(scorevalues) = "numeric"
    rownames(scorevalues) = pairnames
  } else if (LRpair == 'L2R2') {
    scorematrix = LRloop_info$`L2R2 expr_score`
    pairnames = sprintf("%s_%s", scorematrix[,'L2'], scorematrix[,'R2'])
    scorevalues = scorematrix[,c(3:(ncol(scorematrix)-2))]
    class(scorevalues) = "numeric"
    rownames(scorevalues) = pairnames
  }
  
  return(scorevalues)
}