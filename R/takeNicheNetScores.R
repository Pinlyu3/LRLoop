#' Put ligand/receptor_activity scores for expressed ligand-receptor pairs from NicheNet in a specific format that can be further used in the calculation of LRLoop corrected scores.
#' @name takeNicheNetScores
#' @param LRpairs A matrix of LR pairs of interest with columns "from" and "to", in which each row is an LR pair with ligand in the column "from" and receptor in the column "to".
#' @param Basics Result from the function PrepareBasics
#' @param from_to Specify the LRpairs are ligands from celltype1(ct1) to receptors in celltype2(ct2) or ligands from celltype2(ct2) to receptors in celltype1(ct1).
#' @param nichenetscore NicheNet score options, "auroc", "aupr" or "pearson".  Default is "pearson".
#' @param method Options are "Ligand", "Receptor", "LRmean_arithmetic" or "LRmean_geometric".  Default is "Ligand". If method = "Ligand", for each LRpair, only use the ligand activity score. If method = "Receptor", for each LRpair, only use the receptor activity score. If method = "LRmean_arithmetic", for each LRpair, use the arithmetic average of the ligand activity score and the receptor activity score. If method = "LRmean_geometric", for each LRpair, use the geometric average of the ligand activity score and the receptor activity score.
#' @return
#'  A matrix with three columns "from", "to" and "score" with the ligand, receptor and the score value in each column, respectively.
#' @import tidyverse
#' @export



takeNicheNetScores <- function(LRpairs, from_to, Basics, nichenetscore = 'pearson', method = 'Ligand') {
  
  if (from_to == 'ct1_to_ct2') {
    lrexpr = Basics$lr_expr_ct1_to_ct2$bind
    Lactivity = Basics$ligand_activities_matrix_ct1_to_ct2
    Ractivity = Basics$receptor_activities_matrix_ct1_to_ct2
  } else if (from_to == 'ct2_to_ct1') {
    lrexpr = Basics$lr_expr_ct2_to_ct1$bind
    Lactivity = Basics$ligand_activities_matrix_ct2_to_ct1
    Ractivity = Basics$receptor_activities_matrix_ct2_to_ct1
  }
  
  scorematrix = matrix(0, nrow = nrow(LRpairs), ncol = 3)
  colnames(scorematrix) = c("from", "to", "score")
  for (i in 1:nrow(LRpairs)) {
    ligand = LRpairs[i,'from']
    receptor = LRpairs[i,'to']
    idx = which(lrexpr[,'from']==ligand & lrexpr[,'to']==receptor)
    if (length(idx)>0) {
      if (ligand %in% Lactivity[,'test_ligand']) {
        L = as.numeric(Lactivity[Lactivity[,'test_ligand']==ligand, nichenetscore])
      } else {
        L = 0
      }
      if (receptor %in% Ractivity[,'test_receptor']) {
        R = as.numeric(Ractivity[Ractivity[,'test_receptor']==receptor, nichenetscore])
      } else {
        R = 0
      }
      if (method == "Ligand") {
        myscore = L
      } else if (method == 'Receptor') {
        myscore = R
      } else if (method == 'LRmean_arithmetic') {
        myscore = 0.5*(L+R)
      } else if (method == 'LRmean_geometric') {
        myscore = sqrt(L*R)
      }
    } else {
      myscore = 0
    }
    scorematrix[i,] = c(ligand, receptor, myscore)
  }
  
  return(scorematrix)
}