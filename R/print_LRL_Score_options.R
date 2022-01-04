#' Print available options for the function "plotLRL" variable "WhichL1R1L2R2Score"
#'
#' @name print_LRL_Score_options
#' @param LRloop_info The list of the LRloop_network info calculated by function "LRL_info_collection"
#' @import tidyverse
#' @export


print_LRL_Score_options <- function(LRloop_info) {
  LRL_scores = cbind(LRloop_info$`LRloop expr_score`, LRloop_info$`LRloop logFC_score`[,c(5:ncol(LRloop_info$`LRloop logFC_score`))])
  allcols = colnames(LRL_scores)[5:ncol(LRL_scores)]
  options = allcols[substring(allcols,1,8)=="L1R1L2R2"]
  print(options)
}