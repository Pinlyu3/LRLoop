#' Construct a matrix containing receptor-target regulatory probability scores
#'
#' @name construct_receptor_target_matrix
#' @param weighted_networks A list of two weighted networks sig and gr (signaling network and gene regulatory network), in data frame/tibble with columns "from", "to" and "weight".
#' @param receptors A list of receptor gene symbols.
#' @param rtf_cutoff Receptor-tf scores lower than the "rtf_cutoff" quantile will be set to 0.
#' @param damping_factor A number between 0 and 1, the probability that the random walker in the PPR algorithm will continue the walk on the graph.
#' @param secondary_targets "TRUE" or "FALSE", indicates whether putative secondary targets should be included.
#' @param receptors_as_cols "TRUE" or "FALSE", indicates whether receptors should be in columns of the matrix and targets in rows or vice versa.
#' @param remove_direct_links Indicates whether direct ligand-target and receptor-target links in the gene regulatory network should be kept or not. "no": keep links; "ligand": remove direct ligand-target links; "ligand-receptor": remove both direct ligand-target and receptor-target links. 
#' @return
#'  a matrix containing receptor-target scores
#' @import nichenetr tidyverse Seurat dplyr
#' @export


construct_receptor_target_matrix = function(weighted_networks, receptors, rtf_cutoff = 0.99,  
                                          damping_factor = 0.5, secondary_targets = FALSE, receptors_as_cols = TRUE, remove_direct_links = "no") {
  
  requireNamespace("dplyr")
  
  ## workflow; first: give probability score to genes downstream in signaling path starting from receptor.  
  ## second: multiply this matrix with gene regulatory matrix to get the probabilty scores of downstream target genes
  
  # construct receptor-tf matrix
  # remove direct links if required
  if (remove_direct_links != "no"){
    ligands_ = lr_network$from %>% unique()
    receptors_ = lr_network$to %>% unique()
    if (remove_direct_links == "ligand"){
      weighted_networks$gr = weighted_networks$gr %>% filter((from %in% ligands_) == FALSE)
    } else if (remove_direct_links == "ligand-receptor"){
      weighted_networks$gr = weighted_networks$gr %>% filter((from %in% c(ligands_,receptors_)) == FALSE)
    }
  }
  
  rtf_matrix = construct_receptor_tf_matrix(weighted_networks, receptors, rtf_cutoff, damping_factor)
  
  # preparing the gene regulatory matrix
  grn_matrix = construct_tf_target_matrix_noligand(weighted_networks)
  
  # Multiply receptor-tf matrix with tf-target matrix
  receptor_to_target = (rtf_matrix %*% grn_matrix)
  
  # Secondary targets
  if (secondary_targets == TRUE) {
    rtf_matrix = receptor_to_target
    
    if (rtf_cutoff > 0){
      rtf_matrix_TRUE = apply(rtf_matrix,1,function(x){x <= quantile(x,rtf_cutoff)}) %>% t()
      rtf_matrix[rtf_matrix_TRUE] = 0
    }
    
    receptor_to_target_primary = receptor_to_target
    receptor_to_target_secondary = rtf_matrix %*% grn_matrix
    
    # set 0's to number higher than 0 to avoid Inf when inverting in sum
    receptor_to_target_primary[receptor_to_target_primary == 0] = Inf
    receptor_to_target_primary[is.infinite(receptor_to_target_primary)] = min(receptor_to_target_primary)
    
    receptor_to_target_secondary[receptor_to_target_secondary == 0] = Inf
    receptor_to_target_secondary[is.infinite(receptor_to_target_secondary)] = min(receptor_to_target_secondary)
    
    # inverting in sum to emphasize primary targets more (scores secondary targets matrix tend to be higher )
    receptor_to_target = (receptor_to_target_primary **-1 + receptor_to_target_secondary **-1) ** -1
  }
  
  receptor_to_target = receptor_to_target %>% as.matrix()
  
  if (receptors_as_cols == TRUE){
    receptor_to_target = receptor_to_target %>% t()
  }
  
  return(receptor_to_target)
}
