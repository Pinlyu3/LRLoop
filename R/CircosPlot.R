#' Circos plot 
#' @name Circos_plot
#' @param LRscorematrix A matrix with columns "ligand", "receptor" and the condition names where the condition columns contain the LRscores of each ligand-receptor pairs in each condition
#' @param LRscore_conditions A vector of conditions to be used when calculating the maximum or mean LR expr_score for each LR pair acorss conditions
#' @param LRloop_info The list of collected info on the LRloop_network of interest resulted from the function "LRL_info_collection"
#' @param Ligand 'L1' or 'L2'.
#' @param Receptor 'R1' or 'R2'. Remark: Ligand and Receptor should be set in pairs as either Ligand = "L1", Receptor = "R1", or Ligand = "L2", Receptor = "R2"
#' @param WhichLRscore The LR expr_score to be used in the plot, which is 'LRscore_max' (for each LR pair, the maximum LR expr_score across LRscore_conditions) or 'LRscore_mean' (for each LR pair, the average LR expr_score across LRscore_conditions)
#' @param color_L A vector of colors of the ligand clusters, with the name of each element the name of the corresponding cluster
#' @param color_R A vector of colors of the receptor clusters, with the name of each element the name of the corresponding cluster
#' @param color_LR A vector of colors of the LR clusters, with the name of each element the name of the corresponding cluster
#' @param width_same_ligand_cluster Gap width between ligand nodes in the same cluster
#' @param width_different_ligand_cluster Gap width between different ligand clusters 
#' @param width_ligand_receptor Gap width between the ligand nodes and the receptor nodes 
#' @param width_same_receptor_cluster Gap width between receptor nodes in the same cluster
#' @param width_different_receptor_cluster Gap width between different receptor clusters 
#' @param cplotthresh A number. Only plot L-R edges in the circos plot with LRscore > cplotthresh
#' @param cex Font size
#' @import nichenetr tidyverse Seurat circlize igraph RColorBrewer writexl pheatmap dplyr
#' @export



CircosPlot <- function(LRscorematrix, LRscore_conditions, LRloop_info, Ligand, Receptor, WhichLRscore,
                       color_L, color_R, color_LR,
                       width_same_ligand_cluster, width_different_ligand_cluster, width_ligand_receptor, width_same_receptor_cluster, width_different_receptor_cluster,
                       cplotthresh, cex) {
  
  ### Collect LRscores
  LRscore_val = LRscorematrix[,3:ncol(LRscorematrix)]
  mode(LRscore_val) = 'numeric'
  LRscore_val_max = vector()
  LRscore_val_mean = vector()
  for (i in 1:nrow(LRscore_val)) {
    LRscore_val_max[i] = max(LRscore_val[i,LRscore_conditions])
    LRscore_val_mean[i] = mean(LRscore_val[i,LRscore_conditions])
  }
  LRscore = cbind(LRscorematrix, LRscore_val_max, LRscore_val_mean,
                  sprintf("%s_%s", LRscorematrix[,'from'], LRscorematrix[,'to']))
  colnames(LRscore) = c(colnames(LRscorematrix), 'LRscore_max', 'LRscore_mean', 'L_R')
  
  ### Add LRscore to the LRloop network scores
  LRL_scores = LRloop_info$`LRloop expr_score`
  LRL_LR_scores = vector()
  for (i in 1:nrow(LRL_scores)) {
    idx = which(LRscore[,'from']==LRL_scores[i,Ligand] & LRscore[,'to']==LRL_scores[i,Receptor])
    LRL_LR_scores = rbind(LRL_LR_scores, c(LRL_scores[i,], LRscore[idx,3:ncol(LRscore)]))
  }
  colnames(LRL_LR_scores) = c(colnames(LRL_scores), colnames(LRscore)[3:ncol(LRscore)])
  
  LRL_LR_scores_val = LRL_LR_scores[,5:(ncol(LRL_LR_scores)-1)]
  mode(LRL_LR_scores_val) = 'numeric'
  LRLLR_scores = cbind.data.frame(LRL_LR_scores[,c('L1', 'R1', 'L2', 'R2', 'L_R')], LRL_LR_scores_val)
  colnames(LRLLR_scores) = c(c('L1', 'R1', 'L2', 'R2', 'L_R'), colnames(LRL_LR_scores_val))
  
  ######################## Circos data
  LR_df = unique(LRLLR_scores[,c(Ligand, Receptor, WhichLRscore, 'L_R')])
  colnames(LR_df)[1:2] = c('L_name', 'R_name')
  LR_tb = tibble(LR_df)
  
  ### The user pre-defined gene and LR pair clusters
  LRL_clusters_pre = LRloop_info$`L1R1L2R2 clusters`
  LRL_clusters = cbind(LRL_clusters_pre, sprintf("%s_%s", LRL_clusters_pre[,Ligand], LRL_clusters_pre[,Receptor]))
  colnames(LRL_clusters) = c(colnames(LRL_clusters_pre), 'L_R')
  
  # L cluster
  L_clusterset = unique(LRL_clusters[,sprintf("%s_cluster", Ligand)])
  L_cluster_list = list()
  for (i in 1:length(L_clusterset)) {
    L_cluster_list[[i]] = unique(LRL_clusters[LRL_clusters[,sprintf("%s_cluster", Ligand)]==L_clusterset[i],Ligand])
  }
  names(L_cluster_list) = L_clusterset
  L_cluster_vec = vector()
  for (i in 1:length(L_clusterset)) {
    L_cluster_vec = c(L_cluster_vec, rep(L_clusterset[i], length(L_cluster_list[[i]])))
  }
  L_name_vec = vector()
  for (i in 1:length(L_clusterset)) {
    L_name_vec = c(L_name_vec, L_cluster_list[[i]])
  }
  Lcluster_indication_tb = tibble(L_cluster = L_cluster_vec, L_name = L_name_vec)
  
  # R cluster
  R_clusterset = unique(LRL_clusters[,sprintf("%s_cluster", Receptor)])
  R_cluster_list = list()
  for (i in 1:length(R_clusterset)) {
    R_cluster_list[[i]] = unique(LRL_clusters[LRL_clusters[,sprintf("%s_cluster", Receptor)]==R_clusterset[i],Receptor])
  }
  names(R_cluster_list) = R_clusterset
  R_cluster_vec = vector()
  for (i in 1:length(R_clusterset)) {
    R_cluster_vec = c(R_cluster_vec, rep(R_clusterset[i], length(R_cluster_list[[i]])))
  }
  R_name_vec = vector()
  for (i in 1:length(R_clusterset)) {
    R_name_vec = c(R_name_vec, R_cluster_list[[i]])
  }
  Rcluster_indication_tb = tibble(R_cluster = R_cluster_vec, R_name = R_name_vec)
  
  # LR cluster
  LR_clusterset = unique(LRL_clusters[,sprintf("%s%scluster", Ligand, Receptor)])
  LR_cluster_list = list()
  for (i in 1:length(LR_clusterset)) {
    LR_cluster_list[[i]] = unique(LRL_clusters[LRL_clusters[,sprintf("%s%scluster", Ligand, Receptor)]==LR_clusterset[i], 'L_R'])
  }
  names(LR_cluster_list) = LR_clusterset
  LR_cluster_vec = vector()
  for (i in 1:length(LR_clusterset)) {
    LR_cluster_vec = c(LR_cluster_vec, rep(LR_clusterset[i], length(LR_cluster_list[[i]])))
  }
  LR_name_vec = vector()
  for (i in 1:length(LR_clusterset)) {
    LR_name_vec = c(LR_name_vec, LR_cluster_list[[i]])
  }
  LRcluster_indication_tb = tibble(LR_cluster = LR_cluster_vec, L_R = LR_name_vec)
  
  LR_addcluster = unique(inner_join(LR_tb, Lcluster_indication_tb))
  LR_addcluster = unique(inner_join(LR_addcluster, Rcluster_indication_tb))
  LR_addcluster = unique(inner_join(LR_addcluster, LRcluster_indication_tb))
  
  ### Colors
  grid_col_tbl_L = tibble(L_cluster = color_L %>% names(), color_L_cluster = color_L)
  grid_col_tbl_R = tibble(R_cluster = color_R %>% names(), color_R_cluster = color_R)
  grid_col_tbl_LR = tibble(LR_cluster = color_LR %>% names(), color_LR_cluster = color_LR)
  
  LR_addcluster = LR_addcluster %>% mutate(L_name = paste(L_name," ")) 
  
  LR_addcluster = LR_addcluster %>% 
    inner_join(grid_col_tbl_L) %>% 
    inner_join(grid_col_tbl_R) %>% 
    inner_join(grid_col_tbl_LR)
  LR_circle = LR_addcluster %>% select(L_name, R_name, WhichLRscore)
  
  L_color = LR_addcluster %>% distinct(L_name,color_L_cluster)
  grid_L_color = L_color$color_L_cluster %>% set_names(L_color$L_name)
  
  R_color = LR_addcluster %>% distinct(R_name,color_R_cluster)
  grid_R_color = R_color$color_R_cluster %>% set_names(R_color$R_name)
  
  LR_color = LR_addcluster %>% distinct(L_R,color_LR_cluster)
  grid_LR_color = LR_color$color_LR_cluster %>% set_names(LR_color$L_R)
  
  grid_col =c(grid_L_color,grid_R_color)
  
  ### Edge transparency
  if (WhichLRscore == 'LRscore_max') {
    transparency = LR_addcluster %>% mutate(LRscore_max =(LRscore_max-min(LRscore_max))/(max(LRscore_max)-min(LRscore_max))) %>% 
      mutate(transparency = (1.001-LRscore_max)*0.995) %>% .$transparency 
  } else if (WhichLRscore == 'LRscore_mean') {
    transparency = LR_addcluster %>% mutate(LRscore_mean =(LRscore_mean-min(LRscore_mean))/(max(LRscore_mean)-min(LRscore_mean))) %>% 
      mutate(transparency = (1.001-LRscore_mean)*0.995) %>% .$transparency 
  }
  
  ### Order the ligands and the receptors
  L_order = L_name_vec %>% c(paste(.," ")) %>% intersect(LR_addcluster$L_name)
  R_order = R_name_vec %>% c(paste(.,"")) %>% intersect(LR_addcluster$R_name)
  order = c(L_order,R_order)
  
  ### Gaps
  gaps_L = vector()
  for (i in 1:length(L_clusterset)) {
    gaps_L = c(gaps_L, 
               rep(width_same_ligand_cluster,
                   times = (LR_addcluster %>% filter(L_cluster == L_clusterset[i]) %>% distinct(L_name) %>% nrow() -1)),
               width_different_ligand_cluster)
  }
  gaps_L = gaps_L[1:(length(gaps_L)-1)]
  
  gaps_R = vector()
  for (i in 1:length(R_clusterset)) {
    gaps_R = c(gaps_R, 
               rep(width_same_receptor_cluster,
                   times = (LR_addcluster %>% filter(R_cluster == R_clusterset[i]) %>% distinct(R_name) %>% nrow() -1)),
               width_different_receptor_cluster)
  }
  gaps_R = gaps_R[1:(length(gaps_R)-1)]
  
  gaps = c(gaps_L, width_ligand_receptor, gaps_R, width_ligand_receptor)
  
  ### Plot
  LR_circle = LR_addcluster %>% select(L_name,R_name, WhichLRscore, LR_cluster, color_LR_cluster)
  if (WhichLRscore == 'LRscore_max') {
    circos.clear()
    circos.par(gap.degree = gaps)
    chordDiagram(LR_circle, directional = 1, order=order, link.sort = TRUE, link.decreasing = FALSE, col = LR_circle$color_LR_cluster,
                 grid.col = grid_col, transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
                 link.arr.type = "big.arrow", link.visible = LR_circle$LRscore_max > cplotthresh, annotationTrack = "grid", 
                 preAllocateTracks = list(track.height = 0.075))
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = cex)
    }, bg.border = NA)
  } else if (WhichLRscore == 'LRscore_mean') {
    circos.clear()
    circos.par(gap.degree = gaps)
    chordDiagram(LR_circle, directional = 1, order=order, link.sort = TRUE, link.decreasing = FALSE, col = LR_circle$color_LR_cluster,
                 grid.col = grid_col, transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
                 link.arr.type = "big.arrow", link.visible = LR_circle$LRscore_mean > cplotthresh, annotationTrack = "grid", 
                 preAllocateTracks = list(track.height = 0.075))
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = cex)
    }, bg.border = NA)
  }
  
}