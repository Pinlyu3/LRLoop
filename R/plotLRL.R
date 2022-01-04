#' Plot the LRloop network 
#' (A_ and B_ in the LR node lables: Distinguish ligand-from-ct1_receptor-in-ct2 pairs from ligand-from-ct2_receptor-in-ct1 pairs)
#' @name plotLRL
#' @param LRloop_info The list of the LRloop_network info calculated by function "LRL_info_collection"
#' @param WhichL1R1L2R2Score Run the function "print_LRL_Score_options" first to check available options 
#' @param nodecolors  Colors of the A_L1_R1 and B_L2_R2 type nodes in a vector of length two
#' @param labelcolor Label text color
#' @param nodesize Node size
#' @param edgecutoff Only plot edges with L1R1L2R2 LRloopScore greater than edgecutoff
#' @param labelsize Node label size
#' @param edgecolor Edge color
#' @param pheatmapcolor Value of pheatmap variable 'color'
#' @param pheatmap_dist Value of pheatmap variables 'clustering_distance_rows' and 'clustering_distance_cols'. e.g. 'correlation', 'euclidean'
#' @import tidyverse Seurat circlize igraph RColorBrewer writexl pheatmap
#' @export


plotLRL <- function(LRloop_info, WhichL1R1L2R2Score, nodecolors, labelcolor, nodesize, edgecutoff, labelsize, edgecolor, pheatmapcolor, pheatmap_dist) {
  
  LRL_scores = cbind(LRloop_info$`LRloop expr_score`, LRloop_info$`LRloop logFC_score`[,c(5:ncol(LRloop_info$`LRloop logFC_score`))], 
                     LRloop_info$`LRloop nichenet_score`[,c(5:ncol(LRloop_info$`LRloop nichenet_score`))])
  myLRL_edges = cbind.data.frame(sprintf("A_%s_%s", LRL_scores[,'L1'], LRL_scores[,'R1']), sprintf("B_%s_%s", LRL_scores[,'L2'], LRL_scores[,'R2']),
                                 as.numeric(LRL_scores[,WhichL1R1L2R2Score]))
  colnames(myLRL_edges) = c('A_L1_R1', 'B_L2_R2', 'L1R1L2R2Score')
  
  myLRL_nodes = cbind(c(unique(myLRL_edges[,'A_L1_R1']), unique(myLRL_edges[,'B_L2_R2'])), 
                      c(rep('A_L1_R1', length(unique(myLRL_edges[,'A_L1_R1']))), rep('B_L2_R2', length(unique(myLRL_edges[,'B_L2_R2'])))))
  colnames(myLRL_nodes) = c('node', 'type')
  myLRL_nodes = as.data.frame(myLRL_nodes)
  
  ### 
  net = graph_from_data_frame(d=myLRL_edges, vertices=myLRL_nodes, directed=F) 
  names(nodecolors) = c('A_L1_R1', 'B_L2_R2')
  V(net)$label.color = labelcolor
  V(net)$color = nodecolors[V(net)$type]
  V(net)$size = nodesize
  E(net)$width = 1+E(net)$L1R1L2R2Score*1.2
  cut.off = edgecutoff 
  net.sp = delete_edges(net, E(net)[L1R1L2R2Score<cut.off])
  l = layout_with_fr(net)
  plot(net.sp, vertex.label=V(net)$node, vertex.label.cex=labelsize, edge.color=edgecolor, layout=l)
  
  ### Heatmap
  netm = get.adjacency(net, attr="L1R1L2R2Score", sparse=F)
  netm = netm[unique(myLRL_edges[,'A_L1_R1']), unique(myLRL_edges[,'B_L2_R2'])]
  #palf = colorRampPalette(heatmapcolor) 
  #heatmap(netm, Rowv = NA, Colv = NA, col = palf(100), scale="none", margins=c(10,10) )
  
  pheatmap(netm, scale = 'none', cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = T, show_colnames = T, treeheight_col = 0, border_color = FALSE,
           color = pheatmapcolor, clustering_method = 'average', clustering_distance_rows = pheatmap_dist, treeheight_row = 0, legend = FALSE, 
           clustering_distance_cols = pheatmap_dist)
  
  pheatmap(netm, scale = 'none', cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = F, show_colnames = F, treeheight_col = 0, border_color = FALSE,
           color = pheatmapcolor, clustering_method = 'average', clustering_distance_rows = pheatmap_dist, treeheight_row = 0, legend = FALSE,
           clustering_distance_cols = pheatmap_dist)
  
  
}