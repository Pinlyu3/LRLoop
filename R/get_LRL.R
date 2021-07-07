#' Get the LRloop L1_R1-L2_R2 pairs
#' @name get_LRL
#' @param lr_expr_ct2_to_ct1 A ligand-receptor network list resulted from the function "get_expr_lr", signaling from ct2 to ct1
#' @param lr_expr_ct1_to_ct2 A ligand-receptor network list resulted from the function "get_expr_lr", signaling from ct1 to ct2
#' @param ligand_target_matrix_binary_ct2_to_ct1 A matrix of 0s and 1s with target genes in rows and ligands in columns, indicates whether each target gene is identified as a target of each ligand, signaling to/in ct1
#' @param ligand_target_matrix_binary_ct1_to_ct2 A matrix of 0s and 1s with target genes in rows and ligands in columns, indicates whether each target gene is identified as a target of each ligand, signaling to/in ct2
#' @param receptor_target_matrix_binary_ct2_to_ct1 A matrix of 0s and 1s with target genes in rows and receptors in columns, indicates whether each target gene is identified as a target of each receptor, signaling to/in ct1
#' @param receptor_target_matrix_binary_ct1_to_ct2 A matrix of 0s and 1s with target genes in rows and receptors in columns, indicates whether each target gene is identified as a target of each receptor, signaling to/in ct2
#' @return
#'  myLRL: A list of two elements:
#'  myLRL$'L1->L2_L2->L1': LRloop network with L1->L2 and L2->L1 (L2 is a target of L1 and L1 is a target of L2)
#'  myLRL$'R1->L2_R2->L1': LRloop network with R1->L2 and R2->L1 (L2 is a target of R1 and L1 is a target of R2)
#'  myLRL$'eachcondition': LRLoop networks in each condition
#' @import nichenetr tidyverse Seurat
#' @export




get_LRL <- function(lr_expr_ct2_to_ct1, lr_expr_ct1_to_ct2, 
                    ligand_target_matrix_binary_ct2_to_ct1, ligand_target_matrix_binary_ct1_to_ct2,
                    receptor_target_matrix_binary_ct2_to_ct1, receptor_target_matrix_binary_ct1_to_ct2) {
  
  LRL_list = list()
  for (i in 1:length(lr_expr_ct1_to_ct2$eachcondition)) {
    LRL_list[[i]] = make_LRL_list(lr_expr_ct2_to_ct1$eachcondition[[i]], lr_expr_ct1_to_ct2$eachcondition[[i]],
                                      ligand_target_matrix_binary_ct2_to_ct1, ligand_target_matrix_binary_ct1_to_ct2,
                                      receptor_target_matrix_binary_ct2_to_ct1, receptor_target_matrix_binary_ct1_to_ct2)
    LRL_list[[i]] = list(LRL_list[[i]]$`L1->L2 & L2->L1`, LRL_list[[i]]$`R1->L2 & R2->L1`)
    names(LRL_list[[i]]) = c('L1->L2 & L2->L1', 'R1->L2 & R2->L1')
  }
  names(LRL_list) = names(lr_expr_ct1_to_ct2$eachcondition)
  
  LRL_list_bind_L1L2L2L1 = vector()
  LRL_list_bind_R1L2R2L1 = vector()
  for (i in 1:length(LRL_list)) {
    LRL_list_bind_L1L2L2L1 = rbind(LRL_list_bind_L1L2L2L1, LRL_list[[i]]$`L1->L2 & L2->L1`)
    LRL_list_bind_R1L2R2L1 = rbind(LRL_list_bind_R1L2R2L1, LRL_list[[i]]$`R1->L2 & R2->L1`)
  }
  LRL_list_bind_L1L2L2L1 = unique(LRL_list_bind_L1L2L2L1)
  LRL_list_bind_R1L2R2L1 = unique(LRL_list_bind_R1L2R2L1)
  colnames(LRL_list_bind_L1L2L2L1) = c('L2','R2', 'L1', 'R1')
  colnames(LRL_list_bind_R1L2R2L1) = c('L2','R2', 'L1', 'R1')
  
  myLRL = list(LRL_list_bind_L1L2L2L1, LRL_list_bind_R1L2R2L1, LRL_list)
  names(myLRL) = c('L1->L2_L2->L1', 'R1->L2_R2->L1', 'eachcondition')
  
  return(myLRL)
}









