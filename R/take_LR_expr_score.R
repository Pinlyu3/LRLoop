#' Create a dataframe for the expr_scores of L1R1 or L2R2 pairs across all conditions stored in the list resulted from the function 'LRL_info_collection'
#' @name take_LR_expr_score
#' @param LRloop_info A list resulted from the function 'LRL_info_collection' 
#' @param LRpair 'L1R1' or "L2R2'
#' @param LRL_filter "TRUE" or "FALSE.  "TRUE": Take the LRloop filtered version of the LRscores where in each condition, an LRscore is set to 0 if the LR pair does not form any LRloop in that conditon.  "FALSE": Take the LRscores with no LRloop filter applied. 
#' @import tidyverse
#' @export

take_LR_expr_score <- function(LRloop_info, LRpair, LRL_filter = 'FALSE') {
  if (LRpair == 'L1R1' & LRL_filter == "FALSE") {
    scorematrix = LRloop_info$`L1R1 expr_score`
    pairnames = sprintf("%s_%s", scorematrix[,'L1'], scorematrix[,'R1'])
    scorevalues = scorematrix[,c(3:(ncol(scorematrix)-2))]
    class(scorevalues) = "numeric"
    rownames(scorevalues) = pairnames
  } else if (LRpair == 'L2R2' & LRL_filter == "FALSE") {
    scorematrix = LRloop_info$`L2R2 expr_score`
    pairnames = sprintf("%s_%s", scorematrix[,'L2'], scorematrix[,'R2'])
    scorevalues = scorematrix[,c(3:(ncol(scorematrix)-2))]
    class(scorevalues) = "numeric"
    rownames(scorevalues) = pairnames
  } else if (LRpair == 'L1R1' & LRL_filter == "TRUE") {
    scorematrix = LRloop_info$`L1R1 expr_score-LRLoop_filtered`
    pairnames = sprintf("%s_%s", scorematrix[,'L1'], scorematrix[,'R1'])
    scorevalues = scorematrix[,c(3:ncol(scorematrix))]
    class(scorevalues) = "numeric"
    rownames(scorevalues) = pairnames
  } else if (LRpair == 'L2R2' & LRL_filter == "TRUE") {
    scorematrix = LRloop_info$`L2R2 expr_score-LRLoop_filtered`
    pairnames = sprintf("%s_%s", scorematrix[,'L2'], scorematrix[,'R2'])
    scorevalues = scorematrix[,c(3:ncol(scorematrix))]
    class(scorevalues) = "numeric"
    rownames(scorevalues) = pairnames
  }
  
  return(scorevalues)
}