#' Calculate the LRscore of an L-R pair
#'
#' @name LRscore
#' @param L the ligand value used for the calculation of the LRscore
#' @param R the receptor value used for the calculation of the LRscore
#' @param m plays the role of a scalar in the calculation of the LRscore when the calculation method is set to 'scsigr' or individual_scale'
#' @param method The method of calculating the LRscore, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr'. The method 'scsigr' was introduced in SingleCellSignalR.
#' @return
#'  The LRscore of an L-R pair
#' @import tidyverse
#' @export

LRscore = function(L, R, m, method) {
  if (method == 'scsigr') {
    LR = (L*R)^(1/2)/(m+(L*R)^(1/2))
  } else if (method == 'mean') {
    LR = (L+R)/2
  } else if (method == 'individual_scale') {
    LR = (L/(m+L))*(R/(m+R))
  } else if (method == 'individual_scale_exp') {
    LR = exp(-1/L)*exp(-1/R)
  } else if (method == 'product') {
    LR = L*R
  } else if (method =='bias_receptor') {
    LR = L*R/(1-R)
  } else {
    stop("Invalid method")
  }
  
  return(LR)
}





