#' Find the regulatory relationships of all the (L2-R2)-(L1-R1) pairs of interest
#' @name make_LRL_list
#' @param L2R2_network Ligand-receptor network matrix with columns "from" and "to", signaling from ct2 to ct1
#' @param L1R1_network Ligand-receptor network matrix with columns "from" and "to", signaling from ct1 to ct2
#' @param ligand_target_matrix_binary_ct1 A matrix of 0s and 1s with target genes in rows and ligands in columns, indicates whether each target gene is identified as a target of each ligand, signaling to/in ct1
#' @param ligand_target_matrix_binary_ct2 A matrix of 0s and 1s with target genes in rows and ligands in columns, indicates whether each target gene is identified as a target of each ligand, signaling to/in ct2
#' @param receptor_target_matrix_binary_ct1 A matrix of 0s and 1s with target genes in rows and receptors in columns, indicates whether each target gene is identified as a target of each receptor, signaling to/in ct1
#' @param receptor_target_matrix_binary_ct2 A matrix of 0s and 1s with target genes in rows and receptors in columns, indicates whether each target gene is identified as a target of each receptor, signaling to/in ct2
#' @return
#'  p_A_a_list: A list of four elements:
#'  p_A_a[[all_reg_info]]: (L2-R2)-(L1-R1) pairs and the regulatory relationships between these genes (1: TRUE; 0: FALSE).
#'  p_A_a[[L1->L2 & L2->L1]]: (L2-R2)-(L1-R1) pairs that regulate each other, identified based on ligand_target_matrix_binary_ct1 and ligand_target_matrix_binary_ct2  
#'  p_A_a[[R1->L2 & R2->L1]]: (L2-R2)-(L1-R1) pairs that regulate each other, identified based on receptor_target_matrix_binary_ct1 and receptor_target_matrix_binary_ct2  
#'  p_A_a[[L1->L2 & L2->L1 & R1->L2 & R2->L1]]: (L2-R2)-(L1-R1) pairs that regulate each other, identified based on [ligand_target_matrix_binary_ct1 & receptor_target_matrix_binary_ct1] and [ligand_target_matrix_binary_ct2 & receptor_target_matrix_binary_ct2]  
#' @import nichenetr tidyverse
#' @export


make_LRL_list <- function(L2R2_network, L1R1_network,
                            ligand_target_matrix_binary_ct1, ligand_target_matrix_binary_ct2,
                            receptor_target_matrix_binary_ct1, receptor_target_matrix_binary_ct2) {
  
  p_A_a = matrix(0, nrow = nrow(L2R2_network)*nrow(L1R1_network), ncol = 8)
  colnames(p_A_a) = c('L2', 'R2', 'L1', 'R1', 'L1->L2', 'R1->L2', 'L2->L1', 'R2->L1')
  for (i in 1:nrow(L2R2_network)) {
    ligand_a = L2R2_network[i,'from']
    receptor_A = L2R2_network[i,'to']
    for (j in 1:nrow(L1R1_network)) {
      ligand_L = L1R1_network[j,'from']
      receptor_R = L1R1_network[j,'to']
      rowidx = j+(i-1)*nrow(L1R1_network)
      p_A_a[rowidx,'L2'] = ligand_a
      p_A_a[rowidx,'R2'] = receptor_A
      p_A_a[rowidx,'L1'] = ligand_L
      p_A_a[rowidx,'R1'] = receptor_R
      if (ligand_a %in% rownames(ligand_target_matrix_binary_ct2) & ligand_L %in% colnames(ligand_target_matrix_binary_ct2)){
        p_A_a[rowidx,'L1->L2'] = ligand_target_matrix_binary_ct2[ligand_a, ligand_L]
      }
      if (ligand_a %in% rownames(receptor_target_matrix_binary_ct2) & receptor_R %in% colnames(receptor_target_matrix_binary_ct2)) {
        p_A_a[rowidx,'R1->L2'] = receptor_target_matrix_binary_ct2[ligand_a, receptor_R]
      }
      if (ligand_L %in% rownames(ligand_target_matrix_binary_ct1) & ligand_a %in% colnames(ligand_target_matrix_binary_ct1)) {
        p_A_a[rowidx,'L2->L1'] = ligand_target_matrix_binary_ct1[ligand_L, ligand_a]
      }
      if (ligand_L %in% rownames(receptor_target_matrix_binary_ct1) & receptor_A %in% colnames(receptor_target_matrix_binary_ct1)) {
        p_A_a[rowidx,'R2->L1'] = receptor_target_matrix_binary_ct1[ligand_L, receptor_A]
      }
    }
  }
  p_A_a_ligand_target = p_A_a[which(p_A_a[,'L1->L2']=='1' & p_A_a[,'L2->L1']=='1'),1:4]
  p_A_a_receptor_target = p_A_a[which(p_A_a[,'R1->L2']=='1' & p_A_a[,'R2->L1']=='1'),1:4]
  p_A_a_ligand_and_receptor_target = p_A_a[which(p_A_a[,'L1->L2']=='1' & p_A_a[,'L2->L1']=='1' & p_A_a[,'R1->L2']=='1' & p_A_a[,'R2->L1']=='1'),1:4]
  p_A_a_list = list(p_A_a, p_A_a_ligand_target, p_A_a_receptor_target, p_A_a_ligand_and_receptor_target)
  names(p_A_a_list) = c('all_reg_info', 'L1->L2 & L2->L1', 'R1->L2 & R2->L1', 
                        'L1->L2 & L2->L1 & R1->L2 & R2->L1')
  
  return(p_A_a_list)
  
}