#' Get modified LRscores with random LRloop networks
#' @name updateLRscores_rand
#' @param n Number of random LRloop networks to be constructed.
#' @param lr_ct1_to_ct2 A matrix with columns "from" and "to" of ligand-receptor pairs from ct1 to ct2.
#' @param lr_ct2_to_ct1 A matrix with columns "from" and "to" of ligand-receptor pairs from ct2 to ct1.
#' @param LRloops The original LRloop network of all possible LRloops btween ct1 and ct2 of the LRpairs in lr_network.
#' @param LRS Same inputs as those used in updateLRscores.
#' @param LRC Same inputs as those used in updateLRscores.
#' @param from_to Same inputs as those used in updateLRscores.
#' @param LRCform Same inputs as those used in updateLRscores.
#' @param lambda Same inputs as those used in updateLRscores.
#' @param mu Same inputs as those used in updateLRscores.
#' @param k Same inputs as those used in updateLRscores.
#' @return
#'  LRscoresUpdated_nrand: Updated LRscore with n random LRloop networks.
#' @import tidyverse
#' @export


# Get modified LRscores with random LRloop networks
# Inputs:
# n: Number of random LRloop networks to be constructed.
# lr_ct1_to_ct2: A matrix with columns "from" and "to" of ligand-receptor pairs from ct1 to ct2.
# lr_ct2_to_ct1: A matrix with columns "from" and "to" of ligand-receptor pairs from ct2 to ct1.
# LRloops: The original LRloop network of all possible LRloops btween ct1 and ct2 of the LRpairs in lr_network.
# LRS, LRC, from_to, LRCform, lambda, mu, k: Same inputs as those used in updateLRscores.
# Output:
# LRscoresUpdated_nrand: Updated LRscore with n random LRloop networks.

updateLRscores_rand <- function(n, lr_ct1_to_ct2, lr_ct2_to_ct1, LRloops, LRS, LRC, from_to, LRCform, lambda, mu, k) {
  
  LRloops_rand = getnLRloops_rand(n = n, lr_ct1_to_ct2 = lr_ct1_to_ct2, lr_ct2_to_ct1 = lr_ct2_to_ct1, LRloops = LRloops)
  
  LRscoresUpdated_nrand = vector()
  for (i in 1:n) {
    LRscoresUpdated_rand = updateLRscores(LRS = LRS, LRC = LRC, LRloops = LRloops_rand[[i]], from_to = from_to, LRCform = LRCform,
                                            lambda = lambda, mu = mu, k = k)
    LRscoresUpdated_nrand = rbind(LRscoresUpdated_nrand, LRscoresUpdated_rand)
  }
  
  return(LRscoresUpdated_nrand)
}