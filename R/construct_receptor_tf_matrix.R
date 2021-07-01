#' Construct a matrix containing receptor-tf probability scores
#'
#' @name construct_receptor_tf_matrix
#' @param weighted_networks A list of two weighted networks sig and gr (signaling network and gene regulatory network), in data frame/tibble with columns "from", "to" and "weight".
#' @param receptors A list of receptor gene symbols.
#' @param rtf_cutoff Receptor-tf scores lower than the "rtf_cutoff" quantile will be set to 0.
#' @param damping_factor A number between 0 and 1, the probability that the random walker in the PPR algorithm will continue the walk on the graph.
#' @param receptors_as_cols "TRUE" or "FALSE", indicates whether receptors should be in columns of the matrix and targets in rows or vice versa.
#' @return
#'  A matrix containing receptor-tf scores
#' @import nichenetr tidyverse Seurat dplyr
#' @export





# Construct a matrix containing receptor-tf probability scores
# Inputs:
# weighted_networks: A list of two weighted networks sig and gr (signaling network and gene regulatory network), in data frame/tibble with columns "from", "to" and "weight".
# receptors: A list of receptor gene symbols.
# rtf_cutoff: Receptor-tf scores lower than the "rtf_cutoff" quantile will be set to 0. 
# damping_factor: A number between 0 and 1, the probability that the random walker in the PPR algorithm will continue the walk on the graph.
# receptors_as_cols: "TRUE" or "FALSE", indicates whether receptors should be in columns of the matrix and targets in rows or vice versa.
# Output:
# A matrix containing receptor-tf scores

construct_receptor_tf_matrix = function(weighted_networks, receptors, rtf_cutoff = 0.99, damping_factor = 0.5, receptors_as_cols = FALSE) {
  
  requireNamespace("dplyr")
  
  # load in weighted networks
  signaling_network = weighted_networks$sig
  regulatory_network = weighted_networks$gr
  
  # convert ids to numeric for making Matrix::sparseMatrix later on
  allgenes = c(signaling_network$from, signaling_network$to, regulatory_network$from, regulatory_network$to) %>% 
    unique() %>% sort()
  allgenes_integer = allgenes %>% factor() %>% as.numeric()
  allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
  mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
  id2allgenes = mapper(allgenes_id_tbl,"allgenes_integer","allgenes")
  
  signaling_network = signaling_network %>% mutate(from_allgenes = id2allgenes[from], to_allgenes = id2allgenes[to]) %>% 
    arrange(from_allgenes) %>% dplyr::select(from_allgenes,to_allgenes,weight)
  
  # Make Matrix::sparse signaling weighted matrix and graph to apply personalized pagerank
  signaling_network_matrix = Matrix::sparseMatrix(signaling_network$from_allgenes %>% as.integer, 
                                                         signaling_network$to_allgenes %>% as.integer, 
                                                         x=signaling_network$weight %>% as.numeric, 
                                                         dims = c(length(allgenes), length(allgenes)))
  signaling_igraph = igraph::graph_from_adjacency_matrix(signaling_network_matrix, weighted=TRUE, mode="directed")
  # personalized pagerank
  # prepare preference vector E
  E = rep(0,times = length(igraph::V(signaling_igraph)))
  # ppr for every receptor individual
  complete_matrix = lapply(receptors,nichenetr:::PPR_wrapper,E,signaling_igraph,damping_factor,id2allgenes,rtf_cutoff)
  
  rtf_matrix = matrix(unlist(complete_matrix), ncol = length(igraph::V(signaling_igraph)), byrow = TRUE)
  
  rownames(rtf_matrix) = sapply(receptors,function(x){paste0(x,collapse = "-")})
  colnames(rtf_matrix) = allgenes
  
  if (receptors_as_cols == TRUE){
    rtf_matrix = rtf_matrix %>% as.matrix()
    rtf_matrix = rtf_matrix %>% t()
  }
  
  return(rtf_matrix)
}