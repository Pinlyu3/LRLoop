#' Construct a specified number of random LRloop networks between ct1 and ct2.
#' @name getnLRloops_rand
#' @param n Number of random LRloop networks to be constructed.
#' @param lr_ct1_to_ct2 A matrix with columns "from" and "to" of ligand-receptor pairs from ct1 to ct2.
#' @param lr_ct2_to_ct1 A matrix with columns "from" and "to" of ligand-receptor pairs from ct2 to ct1.
#' @param LRloops The original LRloop network of all possible LRloops btween ct1 and ct2 of the LRpairs in lr_network.
#' @return
#'  LRloops_rand: A list of n random LRloop networks.
#' @import tidyverse ggplot2
#' @export


getnLRloops_rand <- function(n, lr_ct1_to_ct2, lr_ct2_to_ct1, LRloops) {
  
  numL1R1 = nrow(unique(LRloops[,c('L1', 'R1')]))
  numL2R2 = nrow(unique(LRloops[,c('L2', 'R2')]))
  LRloops_rand = list()
  for (i in 1:n) {
    L1R1s = lr_ct1_to_ct2[sample(c(1:nrow(lr_ct1_to_ct2)), numL1R1),]
    L2R2s = lr_ct2_to_ct1[sample(c(1:nrow(lr_ct2_to_ct1)), numL2R2),]
    LRloops_rand_L1R1s = L1R1s[sample(x=c(1:nrow(L1R1s)), size = nrow(LRloops), replace = TRUE),]
    LRloops_rand_L2R2s = L2R2s[sample(x=c(1:nrow(L2R2s)), size = nrow(LRloops), replace = TRUE),]
    myLRloops_rand = cbind(LRloops_rand_L2R2s, LRloops_rand_L1R1s)
    colnames(myLRloops_rand) = colnames(LRloops)
    myLRloops_rand = unique(myLRloops_rand)
    LRloops_rand[[i]] = myLRloops_rand
  }
  
  return(LRloops_rand)
}