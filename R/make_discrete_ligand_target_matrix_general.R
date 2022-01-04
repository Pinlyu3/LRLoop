#' Identify the targets of each ligand's cognate receptors (1: TRUE; 0: FALSE)
#' @name make_discrete_ligand_target_matrix_general
#' @param ligands Vector of ligand gene symbols of interest
#' @param receptor_target_matrix_binary Matrix of 0s and 1s with targets in rows and receptors in columns indicating whether each target gene is identified as a target of each receptor.
#' @param lr_network Ligand-receptor network matrix with columns "from" and "to".
#' @return
#'  ligand_target_matrix_general_binary: Matrix of 0s and 1s, indicates whether each gene (in rows) is a target of at least one of the ligand's (in columns) cognate receptors
#' @import nichenetr tidyverse Seurat
#' @export


make_discrete_ligand_target_matrix_general <- function(ligands, receptor_target_matrix_binary, lr_network) {
  
  ligand_target_matrix_general_binary = matrix(0, nrow = nrow(receptor_target_matrix_binary), ncol = length(ligands))
  rownames(ligand_target_matrix_general_binary) = rownames(receptor_target_matrix_binary)
  colnames(ligand_target_matrix_general_binary) = ligands
  for (i in 1:ncol(ligand_target_matrix_general_binary)) {
    ligand = colnames(ligand_target_matrix_general_binary)[i]
    receptors = intersect(lr_network[which(lr_network[,'from']==ligand),'to'], colnames(receptor_target_matrix_binary))
    if (length(receptors)>0) {
      for (j in 1:nrow(ligand_target_matrix_general_binary)) {
        target = rownames(ligand_target_matrix_general_binary)[j]
        ligand_target_matrix_general_binary[j,i] = (sum(receptor_target_matrix_binary[target,receptors])>0)*1
      }
    }
  }
  
  return(ligand_target_matrix_general_binary)
  
}