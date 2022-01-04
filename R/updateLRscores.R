#' Modify the LRscores of the LRs expressed from ct1 to ct2 (from ct2 to ct1) by the LRscores of their LRloop paired LRs from ct2 to ct1 (from ct1 to ct2).
#' @name updateLRscores
#' @param LRS The LRscores to be updated.  Its first two columns 'from' and 'to' are the LRpairs of interest.
#' @param LRC LRScore corrector matrix.  Its first two columns 'from' and 'to' are all the LRpairs in lr_network.
#' @param LRloops A matrix of all possible LRloops (in rows) between ct1 and ct2 of the LRpairs in lr_network.
#' @param from_to Specifies if the LRpairs in LRS are expressed from ct1 to ct2 or from ct2 to ct1.  Options are "ct1_to_ct2" or "ct2_to_ct1".  
#' @param LRCform If an LRpair in LRS has multiple LRloop-partners in LRloops, it specifies how to get a single corrector score for it from the LRscores of these LRloop-partners in LRC. Options are "max", "mean" and "sum".
#' @param lambda Formula parameters in the function L1R1L2R2ScoreHill.
#' @param mu Formula parameters in the function L1R1L2R2ScoreHill.
#' @param k Formula parameters in the function L1R1L2R2ScoreHill.
#' @return
#'  LRSC: LRscore matrix updated from LRS.  Compared to LRS, it has four more columns 'cond_max', 'cond_mean', 'twowaymax' and 'twowaymean'. 
#' @import tidyverse
#' @export


updateLRscores <- function(LRS, LRC, LRloops, from_to, LRCform = 'max', lambda, mu, k) {
  
  LRSC = matrix(0, nrow = nrow(LRS), ncol = ncol(LRS)+4)
  colnames(LRSC) = c(colnames(LRS), 'cond_max', 'cond_mean', 'twowaymax', 'twowaymean')
  conds = colnames(LRS)[3:ncol(LRS)]
  
  for (i in 1:nrow(LRS)) {
    ligand = LRS[i,'from']
    receptor = LRS[i,'to']
    LRSi = as.numeric(LRS[i,conds]) #a vector
    names(LRSi) = conds
    LRSimax = max(LRSi)
    LRSimean = mean(LRSi)
    
    if (from_to == 'ct1_to_ct2') {
      LRloopidx = which(LRloops[,'L1']==ligand & LRloops[,'R1']==receptor)
      if (length(LRloopidx)==0) {
        LRCi = rep(0, length(conds))
      } else {
        LRcorrectorpairs = LRloops[LRloopidx, c('L2', 'R2'), drop=FALSE]
        LRCimatrix = vector()
        for (j in 1:length(LRloopidx)) {
          LRCimatrixj = as.numeric(LRC[LRC[,'from']==LRcorrectorpairs[j,'L2'] & LRC[,'to']==LRcorrectorpairs[j,'R2'],conds,drop=FALSE])
          LRCimatrix = rbind(LRCimatrix, LRCimatrixj)
        }
        if (LRCform == 'max') {
          LRCi = apply(LRCimatrix,2,max)
        } else if (LRCform == 'mean') {
          LRCi = apply(LRCimatrix,2,mean)
        } else if (LRCform == 'sum') {
          LRCi = apply(LRCimatrix,2,sum)
        }
      }
    } else if (from_to == 'ct2_to_ct1') {
      LRloopidx = which(LRloops[,'L2']==ligand & LRloops[,'R2']==receptor)
      if (length(LRloopidx)==0) {
        LRCi = rep(0, length(conds))
      } else {
        LRcorrectorpairs = LRloops[LRloopidx, c('L1', 'R1'), drop=FALSE]
        LRCimatrix = vector()
        for (j in 1:length(LRloopidx)) {
          LRCimatrixj = as.numeric(LRC[LRC[,'from']==LRcorrectorpairs[j,'L1'] & LRC[,'to']==LRcorrectorpairs[j,'R1'],conds,drop=FALSE])
          LRCimatrix = rbind(LRCimatrix, LRCimatrixj)
        }
        if (LRCform == 'max') {
          LRCi = apply(LRCimatrix,2,max)
        } else if (LRCform == 'mean') {
          LRCi = apply(LRCimatrix,2,mean)
        } else if (LRCform == 'sum') {
          LRCi = apply(LRCimatrix,2,sum)
        }
      }
    }
    names(LRCi) = conds
    
    LRCimax = max(LRCi)
    LRCimean = mean(LRCi)
    
    LRSC[i,'from'] = ligand
    LRSC[i,'to'] = receptor
    LRSC[i,'twowaymax'] = L1R1L2R2ScoreHill(myscore = LRSimax, corrector = LRCimax, lambda = lambda, mu = mu, k = k)
    LRSC[i,'twowaymean'] = L1R1L2R2ScoreHill(myscore = LRSimean, corrector = LRCimean, lambda = lambda, mu = mu, k = k)
    for (j in 1:length(conds)) {
      LRSC[i,conds[j]] = L1R1L2R2ScoreHill(myscore = LRSi[conds[j]], corrector = LRCi[conds[j]], lambda = lambda, mu = mu, k = k)
    }
    LRSC[i,'cond_max'] = max(as.numeric(LRSC[i,conds]))
    LRSC[i,'cond_mean'] = mean(as.numeric(LRSC[i,conds]))
  }
  
  return(LRSC)
}