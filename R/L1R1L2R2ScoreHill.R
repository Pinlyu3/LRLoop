#' L1R1L2R2ScoreHill
#' @name L1R1L2R2ScoreHill
#' @param myscore myscore
#' @param corrector corrector
#' @param lambda lambda
#' @param mu mu
#' @param k k
#' @return
#'   z z
#' @import tidyverse ggplot2
#' @export


L1R1L2R2ScoreHill = function(myscore, corrector, lambda, mu, k) {
  x = myscore
  y = corrector
  z = x*(1 + lambda*y^k/(mu^k + y^k))
  return(z)
}