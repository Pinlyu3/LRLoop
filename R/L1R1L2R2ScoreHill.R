#' L1R1L2R2ScoreHill
#' @name L1R1L2R2ScoreHill
#' @param myscore 
#' @param corrector 
#' @param lambda 
#' @param mu 
#' @param k
#' @return
#'   z
#' @import tidyverse ggplot2
#' @export


L1R1L2R2ScoreHill = function(myscore, corrector, lambda, mu, k) {
  x = myscore
  y = corrector
  
  z = x*(1 + lambda*y^k/(mu^k + y^k))
  
  return(z)
}