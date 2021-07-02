#' Calculate the LRscore of an L-R pair
#'
#' @name LRscore
#'
#' @param L A number, the ligand value used for the calculation of the LRscore
#'
#' @param R A number, the receptor value used for the calculation of the LRscore
#'
#' @param m A number, plays the role of a scalar in the calculation of the LRscore when the calculation method is set to 'scsigr' or individual_scale'
#'
#' @param method The method of calculating the LRscore, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr'. The method 'scsigr' was introduced in SingleCellSignalR. Reference: Simon Cabello-Aguilar, Mélissa Alame, Fabien Kon-Sun-Tack, Caroline Fau, Matthieu Lacroix, Jacques Colinge. Nucleic Acids Research, Volume 48, Issue 10, 04 June 2020, Page e55, https://doi.org/10.1093/nar/gkaa183
#'  
#' @return The LRscore of an L-R pair
#'  
#' @import nichenetr tidyverse Seurat ggplot2
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





