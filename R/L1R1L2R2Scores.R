#' Calculate the LRloop score of a L1-R1_L2-R2 pair
#' @name L1R1L2R2Score
#' @param L1 A number, the L1 value used for the calculation of the LRloop score
#' @param R1 A number, the R1 value used for the calculation of the LRloop score
#' @param L2 A number, the L2 value used for the calculation of the LRloop score
#' @param R2 A number, the R2 value used for the calculation of the LRloop score
#' @param m A number, plays the role of a scalar in the calculation of LR scores when LRscore is set to 'scsigr' or individual_scale'
#' @param LRscore The method of calculating the LR scores, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr' 
#' @param LoopScore The method of calculating the LRloop score, available options are "ave_geo" and "ave_alg"
#' @return
#'  The LRloop score of a L1-R1_L2-R2 pair
#' @import nichenetr tidyverse Seurat
#' @export



L1R1L2R2Score = function(L1, R1, L2, R2, m, LRscore, LoopScore) {
  if (LRscore == 'scsigr') {
    L1R1 = (L1*R1)^(1/2)/(m+(L1*R1)^(1/2))
    L2R2 = (L2*R2)^(1/2)/(m+(L2*R2)^(1/2))
  } else if (LRscore == 'mean') {
    L1R1 = (L1+R1)/2
    L2R2 = (L2+R2)/2
  } else if (LRscore == 'individual_scale') {
    L1R1 = (L1/(m+L1))*(R1/(m+R1))
    L2R2 = (L2/(m+L2))*(R2/(m+R2))
  } else if (LRscore == 'individual_scale_exp') {
    L1R1 = exp(-1/L1)*exp(-1/R1)
    L2R2 = exp(-1/L2)*exp(-1/R2)
  } else if (LRscore == 'product') {
    L1R1 = L1*R1
    L2R2 = L2*R2
  } else if (LRscore =='bias_receptor') {
    L1R1 = L1*R1/(1-R1)
    L2R2 = L2*R2/(1-R2)
  } else {
    stop("Invalid LRscore")
  }
  
  if (LoopScore == 'ave_geo') {
    L1R1_L2R2 = (L1R1*L2R2)^(1/2)
  } else if (LoopScore == 'L1R1') {
    L1R1_L2R2 = L1R1
  } else if (LoopScore == 'ave_alg') {
    L1R1_L2R2 = mean(c(L1R1, L2R2))
  } else {
    stop("Invalid LoopScore")
  }
  
  return(L1R1_L2R2)
}


















