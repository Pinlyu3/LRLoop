#' Downweight the hub nodes in the signaling and gene regulatory networks
#'
#' @name apply_hub_corrections_noligand
#' @param weighted_networks A list of three weighted networks lr, sig and gr (ligand-receptor network, signaling network and gene regulatory network), in data frame/tibble with columns "from", "to" and "weight"
#' @param sig_hub A number between 0 and 1. 0: no correction for hubiness; 1: maximal correction for hubiness.
#' @param gr_hub A number between 0 and 1. 0: no correction for hubiness; 1: maximal correction for hubiness.
#' @return
#'  A list of two hubiness-corrected weighted networks sig and gr (signaling network and gene regulatory network), in data frame/tibble with columns "from", "to" and "weight".
#' @import nichenetr tidyverse Seurat
#' @export


apply_hub_corrections_noligand = function(weighted_networks,sig_hub,gr_hub) {
  
  requireNamespace("dplyr")
  
  # load in weighted networks
  signaling_network = weighted_networks$sig
  regulatory_network = weighted_networks$gr
  
  # apply hub correction signaling network
  if (sig_hub > 0){
    signaling_network = signaling_network %>% group_by(to) %>% count(to) %>% ungroup() %>% 
      inner_join(signaling_network, ., by = "to") %>% group_by(from) %>% mutate(weight = weight/(n**sig_hub)) %>% dplyr::select(-n)
  }
  # apply hub correction gene regulatory network
  if (gr_hub > 0){
    regulatory_network = regulatory_network %>% group_by(to) %>% count(to) %>% ungroup() %>% 
      inner_join(regulatory_network, ., by = "to") %>% group_by(from) %>% mutate(weight = weight/(n**gr_hub)) %>% dplyr::select(-n)
  }
  return(list(sig = signaling_network %>% ungroup(), gr = regulatory_network %>% ungroup()))
}