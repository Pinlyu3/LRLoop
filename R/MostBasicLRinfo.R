#' Print: 
#' Number of L1-R1 pairs
#' Number of L1
#' Number of R1
#' Number of L2-R2 pairs
#' List the ligands L1
#' List the receptors R1
#' List the ligands L2
#' List the receptors R2
#' @name MostBasicLRinfo
#' @param LRloop_info The list of the LRloop_network info calculated by function "LRL_info_collection"
#' @import nichenetr tidyverse Seurat
#' @export




MostBasicLRinfo <- function(LRloop_info) {
  print("Number of L1-R1 pairs:")
  print(nrow(unique(LRloop_info$`LRloop expression`[,c('L1', 'R1')])))
  print("Number of L1:")
  print(length(unique(LRloop_info$`LRloop expression`[,'L1'])))
  print("Number of R1:")
  print(length(unique(LRloop_info$`LRloop expression`[,'R1'])))
  print("Number of L2-R2 pairs:")
  print(nrow(unique(LRloop_info$`LRloop expression`[,c('L2', 'R2')])))
  print("Ligands L1:")
  print(unique(LRloop_info$`LRloop expression`[,'L1']))
  print("Receptors R1:")
  print(unique(LRloop_info$`LRloop expression`[,'R1']))
  print("Ligands L2:")
  print(unique(LRloop_info$`LRloop expression`[,'L2']))
  print("Receptors R2:")
  print(unique(LRloop_info$`LRloop expression`[,'R2']))
}