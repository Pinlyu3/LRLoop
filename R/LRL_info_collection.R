#' Collect data on LRloop network relevant ligands and receptors
#'
#' @name LRL_info_collection
#' @param LRloop_network An LRloop network matrix with columns "L1", "R1", "L2" and "R2".
#' @param ave_expr_ct1  Matrix of average expression value of each gene in each condition in ct1/ct2 cells with genes in rows and conditions in columns.
#' @param ave_expr_ct2  Matrix of average expression value of each gene in each condition in ct1/ct2 cells with genes in rows and conditions in columns.
#' @param pct_expr_ct1 Matrix of fraction of ct1/ct2 cells that expressed each gene in each condition with genes in rows and conditions in columns.
#' @param pct_expr_ct2 Matrix of fraction of ct1/ct2 cells that expressed each gene in each condition with genes in rows and conditions in columns.
#' @param DEGinfo_ct1  List of two elements.  The first element DEGinfo_ct1/ct2$DEG is a list where each element is a matrix with genes in rows and at least two colomns "ave_log2FC" and "p_val_adj";  The second element DEGinfo_ct1/ct2$DEgenes is a vector that contains the gene symbols of DEGs in ct1/ct2.
#' @param DEGinfo_ct2  List of two elements.  The first element DEGinfo_ct1/ct2$DEG is a list where each element is a matrix with genes in rows and at least two colomns "ave_log2FC" and "p_val_adj";  The second element DEGinfo_ct1/ct2$DEgenes is a vector that contains the gene symbols of DEGs in ct1/ct2.
#' @param genes_cluster_ct1 User defined ct1/ct2 gene clustering vector (integers 0,1,2,...) with the name of each element the corresponding gene symbol 
#' @param genes_cluster_ct2 User defined ct1/ct2 gene clustering vector (integers 0,1,2,...) with the name of each element the corresponding gene symbol 
#' @param L1R1_cluster User defined ct1 to ct2/ct2 to ct1 ligand-receptor pair clustering vector (integers 0,1,2,...) with the name of each element the corresponding ligand-receptor gene symbols (in the form of L-R)
#' @param L1R2_cluster User defined ct1 to ct2/ct2 to ct1 ligand-receptor pair clustering vector (integers 0,1,2,...) with the name of each element the corresponding ligand-receptor gene symbols (in the form of L-R)
#' @param valuse_ct1 Matrix of the values to be used in the calculation of LRscores and L1R1L2R2 LoopScores with genes in rows and conditions in columns corresponding to ct1/ct2.  
#' @param valuse_ct2 Matrix of the values to be used in the calculation of LRscores and L1R1L2R2 LoopScores with genes in rows and conditions in columns corresponding to ct1/ct2.  
#' @param scalar A number, plays the role of a scalar in the calculation of LRscores when LRscore_method is set to 'scsigr' or individual_scale'.
#' @param ScoreConditionIdx Column indices of valuse_ct1 (same as that of valuse_ct2) to be used in the calculation of the maximum or mean LRscores and LoopScores across conditions
#' @param LRscore_method Method used in the calculation of LRscores, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr' 
#' @param LoopScore_method  Method used in the calculation of L1R1L2R2 LoopScores, available options are "ave_geo" and "ave_alg"
#' @param LRlfcscore_method Method used in the calculation of logFC-based LRscores, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr' 
#' @param LooplfcScore_method Method used in the calculation of logFC-based L1R1L2R2 LoopScores, available options are "ave_geo" and "ave_alg" 
#' @param ligand_activities_matrix_ct1_to_ct2 Ligand/receptor_activities_matrix_ct1_to_ct2 calculated by nichenetr's algorithm
#' @param receptor_activities_matrix_ct1_to_ct2 Ligand/receptor_activities_matrix_ct1_to_ct2 calculated by nichenetr's algorithm
#' @param ligand_activities_matrix_ct2_to_ct1 Ligand/receptor_activities_matrix_ct2_to_ct1 calculated by nichenetr's algorithm
#' @param receptor_activities_matrix_ct2_to_ct1 Ligand/receptor_activities_matrix_ct2_to_ct1 calculated by nichenetr's algorithm
#'
#' @return
#'  LRL_info_collection: a list of matrices:
#'  LRL_info_collection$'LRloop expression': The LRloop network, together with the average expression value of each gene and the fraction of cells expressing each gene in the corresponding cell type and each condition
#'  LRL_info_collection$'L1R1 expression': Identdified ligand-receptor pairs from ct1 to ct2, together with the average expression value of each gene and the fraction of cells expressing each gene in the corresponding cell type and each condition
#'  LRL_info_collection$'L2R2 expression': Identified ligand-receptor pairs from ct2 to ct1, together with the average expression value of each gene and the fraction of cells expressing each gene in the corresponding cell type and each condition
#'  LRL_info_collection$'LRloopDEG': Differential expresssion information of the genes in LRloop_network 
#'  LRL_info_collection$'L1R1DEG': Differential expression information of L1R1 pairs
#'  LRL_info_collection$'L2R2DEG': Differential expression information of L2R2 pairs
#'  LRL_info_collection$'L1R1L2R2 clusters': User defined clusters of L1, R1, L2, R2, L1-R1 pair and L2-R2 pair
#'  LRL_info_collection$'LRloop expr_score': The L1R1L2R2 LoopScores calculated based on valuse_ct1 and valuse_ct2
#'  LRL_info_collection$'LRloop logFC_score': The L1R1L2R2 LoopScores calculated based on the absolute values of the logFC of L1, R1, L2 and R2
#'  LRL_info_collection$'LRloop nichenet_score': The nichenet algorithm based scores of L1, R1, L2 and L2 in each L1R1L2R2 loop
#'  LRL_info_collection$'L1R1 expr_score': The LRscores of L1R1 pairs calculated based on valuse_ct1 and valuse_ct2
#'  LRL_info_collection$'L1R1 logFC_score': The LRscores of L1R1 pairs calculated based on the absolute values of the logFC of L1 and R1
#'  LRL_info_collection$'L1R1 nichenet_score': The nichenet algorithm based scores of L1 and R1 in each L1R1 pair
#'  LRL_info_collection$'L2R2 expr_score': The LRscores of the L2R2 pairs calculated based on valuse_ct2 and valuse_ct1
#'  LRL_info_collection$'L2R2 logFC_score': The LRscores of the L2R2 pairs calculated based on the absolute values of the logFC of L2 and R2
#'  LRL_info_collection$'L2R2 nichenet_score': The nichenet algorithm based scores of L2 and R2 in each L2R2 pair
#'
#' @import nichenetr tidyverse Seurat
#' @export



LRL_info_collection <- function(valuse_ct1, valuse_ct2, scalar, ScoreConditionIdx, LRscore_method, LoopScore_method, LRlfcscore_method, LooplfcScore_method,
                                LRloop_network, 
                                ave_expr_ct1, ave_expr_ct2,
                                pct_expr_ct1, pct_expr_ct2,
                                DEGinfo_ct1, DEGinfo_ct2,
                                genes_cluster_ct1, genes_cluster_ct2, L1R1_cluster, L2R2_cluster,
                                ligand_activities_matrix_ct1_to_ct2, receptor_activities_matrix_ct1_to_ct2,
                                ligand_activities_matrix_ct2_to_ct1, receptor_activities_matrix_ct2_to_ct1) {
  
  LRloop_info = list()
  
  ################################################################## 'LRloop expression' ##################################################################
  LRloop_info[[1]] = matrix(0, nrow = nrow(LRloop_network), ncol = (2*ncol(ave_expr_ct1)+2*ncol(ave_expr_ct2)+2*ncol(pct_expr_ct1)+2*ncol(pct_expr_ct2)+4))
  for (i in 1:nrow(LRloop_network)) {
    L1 = LRloop_network[i,'L1']
    R1 = LRloop_network[i,'R1']
    L2 = LRloop_network[i,'L2']
    R2 = LRloop_network[i,'R2']
    L1_ave_expr_vec = ave_expr_ct1[L1,]
    R1_ave_expr_vec = ave_expr_ct2[R1,]
    L2_ave_expr_vec = ave_expr_ct2[L2,]
    R2_ave_expr_vec = ave_expr_ct1[R2,]
    L1_pct_expr_vec = pct_expr_ct1[L1,]
    R1_pct_expr_vec = pct_expr_ct2[R1,]
    L2_pct_expr_vec = pct_expr_ct2[L2,]
    R2_pct_expr_vec = pct_expr_ct1[R2,]
    LRloop_info[[1]][i,] = c(c(L1,R1,L2,R2,
                               L1_ave_expr_vec,R1_ave_expr_vec,L2_ave_expr_vec,R2_ave_expr_vec,
                               L1_pct_expr_vec,R1_pct_expr_vec,L2_pct_expr_vec,R2_pct_expr_vec))
  }
  colnames(LRloop_info[[1]]) = c("L1", "R1", "L2", "R2", 
                                 sprintf("L1_ave_expr_%s", colnames(ave_expr_ct1)), sprintf("R1_ave_expr_%s", colnames(ave_expr_ct2)),
                                 sprintf("L2_ave_expr_%s", colnames(ave_expr_ct2)), sprintf("R2_ave_expr_%s", colnames(ave_expr_ct1)),
                                 sprintf("L1_pct_expr_%s", colnames(pct_expr_ct1)), sprintf("R1_pct_expr_%s", colnames(pct_expr_ct2)),
                                 sprintf("L2_pct_expr_%s", colnames(pct_expr_ct2)), sprintf("R2_pct_expr_%s", colnames(pct_expr_ct1)))
  
  ################################################################## 'L1R1 expression' ##################################################################
  LRloop_info[[2]] = unique(LRloop_info[[1]][,c("L1", "R1", 
                                                sprintf("L1_ave_expr_%s", colnames(ave_expr_ct1)), sprintf("R1_ave_expr_%s", colnames(ave_expr_ct2)),
                                                sprintf("L1_pct_expr_%s", colnames(pct_expr_ct1)), sprintf("R1_pct_expr_%s", colnames(pct_expr_ct2)))])
  
  ################################################################## 'L2R2 expression' ##################################################################
  LRloop_info[[3]] = unique(LRloop_info[[1]][,c("L2", "R2", 
                                                sprintf("L2_ave_expr_%s", colnames(ave_expr_ct2)), sprintf("R2_ave_expr_%s", colnames(ave_expr_ct1)),
                                                sprintf("L2_pct_expr_%s", colnames(pct_expr_ct2)), sprintf("R2_pct_expr_%s", colnames(pct_expr_ct1)))])
  
  ################################################################## 'LRloopDEG' ##################################################################
  LRloop_info[[4]] = matrix(0, nrow = nrow(LRloop_network), ncol = (8+
                                                                      2*length(DEGinfo_ct1$DEG)*ncol(DEGinfo_ct1$DEG[[1]])+
                                                                      2*length(DEGinfo_ct2$DEG)*ncol(DEGinfo_ct2$DEG[[1]])))
  for (i in 1:nrow(LRloop_network)) {
    L1 = LRloop_network[i,'L1']
    R1 = LRloop_network[i,'R1']
    L2 = LRloop_network[i,'L2']
    R2 = LRloop_network[i,'R2']
    L1_diff = (L1 %in% DEGinfo_ct1$DEgenes)
    R1_diff = (R1 %in% DEGinfo_ct2$DEgenes)
    L2_diff = (L2 %in% DEGinfo_ct2$DEgenes)
    R2_diff = (R2 %in% DEGinfo_ct1$DEgenes)
    L1_vec = vector()
    R2_vec = vector()
    for (j in 1:length(DEGinfo_ct1$DEG)) {
      L1_vec = c(L1_vec, as.matrix(DEGinfo_ct1$DEG[[j]][L1,]))          
      R2_vec = c(R2_vec, as.matrix(DEGinfo_ct1$DEG[[j]][R2,]))
    }
    R1_vec = vector()
    L2_vec = vector()
    for (j in 1:length(DEGinfo_ct2$DEG)) {
      R1_vec = c(R1_vec, as.matrix(DEGinfo_ct2$DEG[[j]][R1,]))          
      L2_vec = c(L2_vec, as.matrix(DEGinfo_ct2$DEG[[j]][L2,]))
    }
    LRloop_info[[4]][i,] = c(L1, R1, L2, R2,
                             L1_diff, R1_diff, L2_diff, R2_diff,
                             L1_vec, R1_vec, L2_vec, R2_vec)
  }
  diffct1names = vector()
  for (j in 1:length(DEGinfo_ct1$DEG)) {
    diffct1names = c(diffct1names, sprintf("%s_%s", colnames(DEGinfo_ct1$DEG[[1]]), names(DEGinfo_ct1$DEG)[j]))
  }
  diffct2names = vector()
  for (j in 1:length(DEGinfo_ct2$DEG)) {
    diffct2names = c(diffct2names, sprintf("%s_%s", colnames(DEGinfo_ct2$DEG[[1]]), names(DEGinfo_ct2$DEG)[j]))
  }
  colnames(LRloop_info[[4]]) = c("L1", "R1", "L2", "R2",
                                 "L1_diff", "R1_diff", "L2_diff", "R2_diff",
                                 sprintf("L1_%s", diffct1names), sprintf("R1_%s", diffct2names), 
                                 sprintf("L2_%s", diffct2names), sprintf("R2_%s", diffct1names))
  
  ################################################################## 'L1R1DEG' ##################################################################
  LRloop_info[[5]] = unique(LRloop_info[[4]][,c("L1", "R1",
                                                "L1_diff", "R1_diff",
                                                sprintf("L1_%s", diffct1names), sprintf("R1_%s", diffct2names))])
  
  ################################################################## 'L2R2DEG' ##################################################################
  LRloop_info[[6]] = unique(LRloop_info[[4]][,c("L2", "R2",
                                                "L2_diff", "R2_diff",
                                                sprintf("L2_%s", diffct2names), sprintf("R2_%s", diffct1names))])
  
  ################################################################## 'L1R1L2R2 clusters' ##################################################################
  LRloop_info[[7]] = matrix(0, nrow = nrow(LRloop_network), ncol = 10)
  for (i in 1:nrow(LRloop_network)) {
    L1 = LRloop_network[i,'L1']
    R1 = LRloop_network[i,'R1']
    L2 = LRloop_network[i,'L2']
    R2 = LRloop_network[i,'R2']
    L1R1 = sprintf("%s_%s", L1, R1)
    L2R2 = sprintf("%s_%s", L2, R2)
    L1_cluster = genes_cluster_ct1[L1]
    R1_cluster = genes_cluster_ct2[R1]
    L2_cluster = genes_cluster_ct2[L2]
    R2_cluster = genes_cluster_ct1[R2]
    L1R1cluster = L1R1_cluster[L1R1]
    L2R2cluster = L2R2_cluster[L2R2]
    LRloop_info[[7]][i,] = c(L1,R1,L2,R2,L1_cluster, R1_cluster, L2_cluster, R2_cluster, L1R1cluster, L2R2cluster)
  }
  colnames(LRloop_info[[7]]) = c("L1", "R1", "L2", "R2", 'L1_cluster', 'R1_cluster', 'L2_cluster', 'R2_cluster', 'L1R1cluster', 'L2R2cluster')
  
  ################################################################## 'LRloop expr_score' ##################################################################
  LRloop_info[[8]] = matrix(0, nrow = nrow(LRloop_network), ncol = (4+3*ncol(valuse_ct1)+6))
  for (i in 1:nrow(LRloop_network)) {
    L1 = LRloop_network[i,'L1']
    R1 = LRloop_network[i,'R1']
    L2 = LRloop_network[i,'L2']
    R2 = LRloop_network[i,'R2']
    
    L1R1L2R2_expScore = vector()
    L1R1_expscore = vector()
    L2R2_expscore = vector()
    for (j in 1:ncol(valuse_ct1)) {
      cond = colnames(valuse_ct1)[j]
      L1R1L2R2_expScore[j] = L1R1L2R2Score(valuse_ct1[L1,cond], valuse_ct2[R1,cond], valuse_ct2[L2,cond], valuse_ct1[R2,cond],
                                           scalar, LRscore_method, LoopScore_method)
      L1R1_expscore[j] = LRscore(valuse_ct1[L1,cond], valuse_ct2[R1,cond], scalar, LRscore_method)
      L2R2_expscore[j] = LRscore(valuse_ct2[L2,cond], valuse_ct1[R2,cond], scalar, LRscore_method)
    }
    L1R1L2R2_expScore_max = max(L1R1L2R2_expScore[ScoreConditionIdx])
    L1R1L2R2_expScore_mean = mean(L1R1L2R2_expScore[ScoreConditionIdx])
    L1R1_expscore_max = max(L1R1_expscore[ScoreConditionIdx])
    L1R1_expscore_mean = mean(L1R1_expscore[ScoreConditionIdx])
    L2R2_expscore_max = max(L2R2_expscore[ScoreConditionIdx])
    L2R2_expscore_mean = mean(L2R2_expscore[ScoreConditionIdx])
    
    LRloop_info[[8]][i,] = c(L1, R1, L2, R2,
                             L1R1L2R2_expScore, L1R1L2R2_expScore_max, L1R1L2R2_expScore_mean,
                             L1R1_expscore, L1R1_expscore_max, L1R1_expscore_mean,
                             L2R2_expscore, L2R2_expscore_max, L2R2_expscore_mean)
  }
  colnames(LRloop_info[[8]]) = c("L1", "R1", "L2", "R2",
                                 sprintf("L1R1L2R2_exprScore_%s", colnames(valuse_ct1)), "L1R1L2R2_exprScore_max", "L1R1L2R2_exprScore_mean",
                                 sprintf("L1R1_exprscore_%s", colnames(valuse_ct1)), "L1R1_exprscore_max", "L1R1_exprscore_mean",
                                 sprintf("L2R2_exprscore_%s", colnames(valuse_ct1)), "L2R2_exprscore_max", "L2R2_exprscore_mean")
  
  ################################################################## 'L1R1 expr_score' ##################################################################
  LRloop_info[[9]] = unique(LRloop_info[[8]][,c("L1", "R1",
                                                sprintf("L1R1_exprscore_%s", colnames(valuse_ct1)), "L1R1_exprscore_max", "L1R1_exprscore_mean")])
  
  ################################################################## 'L2R2 expr_score' ##################################################################
  LRloop_info[[10]] = unique(LRloop_info[[8]][,c("L2", "R2",
                                                 sprintf("L2R2_exprscore_%s", colnames(valuse_ct1)), "L2R2_exprscore_max", "L2R2_exprscore_mean")])
  
  ################################################################## 'LRloop logFC_score' ##################################################################
  LRloop_info[[11]] = matrix(0, nrow = nrow(LRloop_network), ncol = (4+3*length(DEGinfo_ct1$DEG)+12))
  for (i in 1:nrow(LRloop_network)) {
    L1 = LRloop_network[i,'L1']
    R1 = LRloop_network[i,'R1']
    L2 = LRloop_network[i,'L2']
    R2 = LRloop_network[i,'R2']
    
    L1_lfc_vec = vector()
    for (j in 1:length(DEGinfo_ct1$DEG)) {
      if (!is.na(DEGinfo_ct1$DEG[[j]][L1,'avg_log2FC'])) {
        L1_lfc_vec[j] = abs(DEGinfo_ct1$DEG[[j]][L1,'avg_log2FC'])
      } else {
        L1_lfc_vec[j] = 0
      }
    }
    names(L1_lfc_vec) = names(DEGinfo_ct1$DEG)
    
    R1_lfc_vec = vector()
    for (j in 1:length(DEGinfo_ct2$DEG)) {
      if (!is.na(DEGinfo_ct2$DEG[[j]][R1,'avg_log2FC'])) {
        R1_lfc_vec[j] = abs(DEGinfo_ct2$DEG[[j]][R1,'avg_log2FC'])
      } else {
        R1_lfc_vec[j] = 0
      }
    }
    names(R1_lfc_vec) = names(DEGinfo_ct2$DEG)
    
    L2_lfc_vec = vector()
    for (j in 1:length(DEGinfo_ct2$DEG)) {
      if (!is.na(DEGinfo_ct2$DEG[[j]][L2,'avg_log2FC'])) {
        L2_lfc_vec[j] = abs(DEGinfo_ct2$DEG[[j]][L2,'avg_log2FC'])
      } else {
        L2_lfc_vec[j] = 0
      }
    }
    names(L2_lfc_vec) = names(DEGinfo_ct2$DEG)
    
    R2_lfc_vec = vector()
    for (j in 1:length(DEGinfo_ct1$DEG)) {
      if (!is.na(DEGinfo_ct1$DEG[[j]][R2,'avg_log2FC'])) {
        R2_lfc_vec[j] = abs(DEGinfo_ct1$DEG[[j]][R2,'avg_log2FC'])
      } else {
        R2_lfc_vec[j] = 0
      }
    }
    names(R2_lfc_vec) = names(DEGinfo_ct1$DEG)
    
    L1_maxlfc = max(L1_lfc_vec)
    L1_meanlfc = mean(L1_lfc_vec)
    R1_maxlfc = max(R1_lfc_vec)
    R1_meanlfc = mean(R1_lfc_vec)
    L2_maxlfc = max(L2_lfc_vec)
    L2_meanlfc = mean(L2_lfc_vec)
    R2_maxlfc = max(R2_lfc_vec)
    R2_meanlfc = mean(R2_lfc_vec)
    
    L1R1L2R2_lfcScore = vector()
    L1R1_lfcscore = vector()
    L2R2_lfcscore = vector()
    for (j in 1:length(DEGinfo_ct1$DEG)) {
      cond = names(DEGinfo_ct1$DEG)[j]
      L1R1L2R2_lfcScore[j] = L1R1L2R2Score(L1_lfc_vec[cond], R1_lfc_vec[cond], L2_lfc_vec[cond], R2_lfc_vec[cond],
                                           scalar, LRlfcscore_method, LooplfcScore_method)
      L1R1_lfcscore[j] = LRscore(L1_lfc_vec[cond], R1_lfc_vec[cond], scalar, LRlfcscore_method)
      L2R2_lfcscore[j] = LRscore(L2_lfc_vec[cond], R2_lfc_vec[cond], scalar, LRlfcscore_method)
    }
    
    L1R1L2R2_lfcScore_max = max(L1R1L2R2_lfcScore)
    L1R1L2R2_lfcScore_mean = mean(L1R1L2R2_lfcScore)
    L1R1L2R2_maxlfcbasedScore = L1R1L2R2Score(L1_maxlfc, R1_maxlfc, L2_maxlfc, R2_maxlfc,
                                              scalar, LRlfcscore_method, LooplfcScore_method)
    L1R1L2R2_meanlfcbasedScore = L1R1L2R2Score(L1_meanlfc, R1_meanlfc, L2_meanlfc, R2_meanlfc,
                                               scalar, LRlfcscore_method, LooplfcScore_method)
    
    L1R1_lfcscore_max = max(L1R1_lfcscore)
    L1R1_lfcscore_mean = mean(L1R1_lfcscore)
    L1R1_maxlfcbasedscore = LRscore(L1_maxlfc, R1_maxlfc, scalar, LRlfcscore_method)
    L1R1_meanlfcbasedscore = LRscore(L1_meanlfc, R1_meanlfc, scalar, LRlfcscore_method)
    
    L2R2_lfcscore_max = max(L2R2_lfcscore)
    L2R2_lfcscore_mean = mean(L2R2_lfcscore)
    L2R2_maxlfcbasedscore = LRscore(L2_maxlfc, R2_maxlfc, scalar, LRlfcscore_method)
    L2R2_meanlfcbasedscore = LRscore(L2_meanlfc, R2_meanlfc, scalar, LRlfcscore_method)
    
    LRloop_info[[11]][i,] = c(L1, R1, L2, R2,
                              L1R1L2R2_lfcScore, L1R1L2R2_lfcScore_max, L1R1L2R2_lfcScore_mean, L1R1L2R2_maxlfcbasedScore, L1R1L2R2_meanlfcbasedScore,
                              L1R1_lfcscore, L1R1_lfcscore_max, L1R1_lfcscore_mean, L1R1_maxlfcbasedscore, L1R1_meanlfcbasedscore,
                              L2R2_lfcscore, L2R2_lfcscore_max, L2R2_lfcscore_mean, L2R2_maxlfcbasedscore, L2R2_meanlfcbasedscore)
    
  }
  colnames(LRloop_info[[11]]) = c("L1", "R1", "L2", "R2",
                                  sprintf("L1R1L2R2_logFC_Score_%s", names(DEGinfo_ct1$DEG)), 
                                  "L1R1L2R2_logFC_Score_max", "L1R1L2R2_logFC_Score_mean", "L1R1L2R2_maxlogFC_basedScore", "L1R1L2R2_meanlogFC_basedScore",
                                  sprintf("L1R1_logFC_score_%s", names(DEGinfo_ct1$DEG)),
                                  "L1R1_logFC_score_max", "L1R1_logFC_score_mean", "L1R1_maxlogFC_basedscore", "L1R1_meanlogFC_basedscore",
                                  sprintf("L2R2_logFC_score_%s", names(DEGinfo_ct1$DEG)),
                                  "L2R2_logFC_score_max", "L2R2_logFC_score_mean", "L2R2_maxlogFC_basedscore", "L2R2_meanlogFC_basedscore")
  
  ################################################################## 'L1R1 logFC_score' ##################################################################
  LRloop_info[[12]] = unique(LRloop_info[[11]][,c("L1", "R1", 
                                                  sprintf("L1R1_logFC_score_%s", names(DEGinfo_ct1$DEG)),
                                                  "L1R1_logFC_score_max", "L1R1_logFC_score_mean", "L1R1_maxlogFC_basedscore", "L1R1_meanlogFC_basedscore")])
  
  ################################################################## 'L2R2 logFC_score' ##################################################################
  LRloop_info[[13]] = unique(LRloop_info[[11]][,c("L2", "R2", 
                                                  sprintf("L2R2_logFC_score_%s", names(DEGinfo_ct1$DEG)),
                                                  "L2R2_logFC_score_max", "L2R2_logFC_score_mean", "L2R2_maxlogFC_basedscore", "L2R2_meanlogFC_basedscore")])
  
  ################################################################## 'LRloop nichenet_score' ##################################################################
  LRloop_info[[14]] = matrix(0, nrow = nrow(LRloop_network), ncol = 8)
  for (i in 1:nrow(LRloop_network)) {
    L1 = LRloop_network[i,'L1']
    R1 = LRloop_network[i,'R1']
    L2 = LRloop_network[i,'L2']
    R2 = LRloop_network[i,'R2']
    
    if (L1 %in% ligand_activities_matrix_ct1_to_ct2[,'test_ligand']) {
      L1_nichenetscore = ligand_activities_matrix_ct1_to_ct2[ligand_activities_matrix_ct1_to_ct2[,'test_ligand']==L1,'pearson']
    } else {
      L1_nichenetscore = 0
    }
    if (R1 %in% receptor_activities_matrix_ct1_to_ct2[,'test_receptor']) {
      R1_nichenetscore = receptor_activities_matrix_ct1_to_ct2[receptor_activities_matrix_ct1_to_ct2[,'test_receptor']==R1,'pearson']
    } else {
      R1_nichenetscore = 0
    }
    if (L2 %in% ligand_activities_matrix_ct2_to_ct1[,'test_ligand']) {
      L2_nichenetscore = ligand_activities_matrix_ct2_to_ct1[ligand_activities_matrix_ct2_to_ct1[,'test_ligand']==L2,'pearson']
    } else {
      L2_nichenetscore = 0
    }
    if (R2 %in% receptor_activities_matrix_ct2_to_ct1[,'test_receptor']) {
      R2_nichenetscore = receptor_activities_matrix_ct2_to_ct1[receptor_activities_matrix_ct2_to_ct1[,'test_receptor']==R2,'pearson']
    } else {
      R2_nichenetscore = 0
    }
    
    LRloop_info[[14]][i,] = c(L1, R1, L2, R2, L1_nichenetscore, R1_nichenetscore, L2_nichenetscore, R2_nichenetscore)
  }
  colnames(LRloop_info[[14]]) = c("L1", "R1", "L2", "R2", "L1_nichenetscore", "R1_nichenetscore", "L2_nichenetscore", "R2_nichenetscore")
  
  ################################################################## 'L1R1 nichenet_score' ##################################################################
  LRloop_info[[15]] = unique(LRloop_info[[14]][,c("L1", "R1", "L1_nichenetscore", "R1_nichenetscore")])
  
  ################################################################## 'L2R2 nichenet_score' ##################################################################
  LRloop_info[[16]] = unique(LRloop_info[[14]][,c("L2", "R2", "L2_nichenetscore", "R2_nichenetscore")])
  
  #############################################################################################################################################################
  names(LRloop_info) = c("LRloop expression", "L1R1 expression", "L2R2 expression", "LRloopDEG", "L1R1DEG", "L2R2DEG", 
                         "L1R1L2R2 clusters", "LRloop expr_score", "L1R1 expr_score", "L2R2 expr_score", 
                         "LRloop logFC_score", "L1R1 logFC_score", "L2R2 logFC_score",
                         "LRloop nichenet_score", "L1R1 nichenet_score", "L2R2 nichenet_score")
  
  return(LRloop_info)
}





