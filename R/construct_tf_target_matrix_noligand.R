#' Convert gene regulatory network into a matrix.
#'
#' @name construct_tf_target_matrix_noligand
#' @param weighted_networks the Seurat object of interest
#' @param tfs_as_cols vectors of idents
#' @param standalone_output vectors of idents
#'
#' @return
#'  A matrix of tf-target regulatory weights
#' @import nichenetr tidyverse Seurat dplyr
#' @export



# Convert gene regulatory network into a matrix.
# Inputs:
# weighted_networks: A list of two weighted networks sig and gr (signaling network and gene regulatory network), in data frame/tibble with columns "from", "to" and "weight".
# tfs_as_cols: "TRUE" or "FALSE", indicates whether regulators should be in columns of the matrix and targets in rows or vice versa.
# standalone_output: "TRUE" or "FALSE", indicates whether to keep only regulators with gene regulatory interactions.
# Output:
# A matrix of tf-target regulatory weights

construct_tf_target_matrix_noligand = function(weighted_networks, tfs_as_cols = FALSE, standalone_output = FALSE) {
  
  requireNamespace("dplyr")
  
  # load in weighted networks
  signaling_network = weighted_networks$sig
  regulatory_network = weighted_networks$gr
  
  # convert ids to numeric for making Matrix::sparseMatrix later on
  allgenes = c(signaling_network$from, signaling_network$to, regulatory_network$from, regulatory_network$to) %>% unique() %>% sort()
  allgenes_integer = allgenes %>% factor() %>% as.numeric()
  allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
  mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
  id2allgenes = mapper(allgenes_id_tbl,"allgenes_integer","allgenes")
  
  regulatory_network = regulatory_network %>% mutate(from_allgenes = id2allgenes[from], to_allgenes = id2allgenes[to]) %>% 
    arrange(from_allgenes) %>% dplyr::select(from_allgenes,to_allgenes,weight)
  grn_matrix = Matrix::sparseMatrix(regulatory_network$from_allgenes %>% as.integer, 
                                    regulatory_network$to_allgenes %>% as.integer, 
                                    x=regulatory_network$weight %>% as.numeric, 
                                    dims = c(length(allgenes), length(allgenes)))
  
  rownames(grn_matrix) = allgenes
  colnames(grn_matrix) = allgenes
  if (standalone_output == TRUE){
    # keep only regulators with gene regulatory interactions
    regulators = weighted_networks$gr %>% .$from %>% unique()
    grn_matrix = grn_matrix[regulators,]
    
  }
  if (tfs_as_cols == TRUE){
    grn_matrix = grn_matrix %>% as.matrix()
    grn_matrix = grn_matrix %>% t()
  }
  
  return(grn_matrix)
}