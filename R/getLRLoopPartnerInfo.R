#' Get the number and list of LRLoop partners of LRpairs of interest.
#' @name getLRLoopPartnerInfo
#' @param LRpairs A matirx of the LRpairs of interes with two columns "from" and "to".
#' @param from_to "ct1_to_ct2" or "ct2_to_ct1" specifying if the LRpairs are from ct1 to ct2 or from ct2 to ct1.
#' @param LRLoop_network the LRLoop network between ct1 and ct2.
#' @return
#'  LR_LRLoopPartnerInfo: A matrix with four columns "from", "to", "num_LRLoopPartner", "LRLoopPartner", in "from" and "to" are the LRpairs of interest, in "num_LRLoopPartner" and "LRLoopPartner" are the number and list of LRLoop partners of the LRpairs.  When "num_LRLoopPartner" is 0, the corresponding "LRLoopPartner" is "NA".
#' @import tidyverse Seurat
#' @export



getLRLoopPartnerInfo <- function(LRpairs, from_to, LRLoop_network) {
  LR_LRLoopPartnerInfo = matrix(0, nrow = nrow(LRpairs), ncol = 4)
  colnames(LR_LRLoopPartnerInfo) = c('from','to','num_LRLoopPartner','LRLoopPartner')
  if (from_to == 'ct1_to_ct2') {
    for (i in 1:nrow(LRpairs)) {
      ligand = LRpairs[i,'from']
      receptor = LRpairs[i,'to']
      idx = which(LRLoop_network[,'L1']==ligand & LRLoop_network[,'R1']==receptor)
      numLRLoopPartner = length(idx)
      if (length(idx)>0) {
        LRLoopPartner = 'L2R2'
        for (j in 1:length(idx)) {
          LRLoopPartner = sprintf("%s_%s-%s", LRLoopPartner, LRLoop_network[idx[j],'L2'], LRLoop_network[idx[j],'R2'])
        }
      } else {
        LRLoopPartner = 'NA'
      }
      LR_LRLoopPartnerInfo[i,] = c(ligand, receptor, numLRLoopPartner, LRLoopPartner)
    }
  } else if (from_to == 'ct2_to_ct1') {
    for (i in 1:nrow(LRpairs)) {
      ligand = LRpairs[i,'from']
      receptor = LRpairs[i,'to']
      idx = which(LRLoop_network[,'L2']==ligand & LRLoop_network[,'R2']==receptor)
      numLRLoopPartner = length(idx)
      if (length(idx)>0) {
        LRLoopPartner = 'L1R1'
        for (j in 1:length(idx)) {
          LRLoopPartner = sprintf("%s_%s-%s", LRLoopPartner, LRLoop_network[idx[j],'L1'], LRLoop_network[idx[j],'R1'])
        }
      } else {
        LRLoopPartner = 'NA'
      }
      LR_LRLoopPartnerInfo[i,] = c(ligand, receptor, numLRLoopPartner, LRLoopPartner)
    }
  }
  
  return(LR_LRLoopPartnerInfo)
}