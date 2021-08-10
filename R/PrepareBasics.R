#' Get a list of some basic data info
#' @name PrepareBasics
#' @param ct1obj  Seurat objects of cell type1 and cell type2
#' @param ct2obj  Seurat objects of cell type1 and cell type2
#' @param min_pct Cutoff value on the detection rate of genes in each cell type, each condition for the identification of "expressed" genes. A number between 0 and 1. For example, min_pct = 0.1 means that for each cell type and condition, the genes that are detected in at least 10 percent of the cells are considered "expressed".
#' @param geneset_ct1 Vectors of gene symbols considered as the target gene set of interest (for example, the differentailly expressed genes) in cell type1/cell type2 for the calculation of the nichenetr based ligand/receptor activity scores.  
#' @param geneset_ct2 Vectors of gene symbols considered as the target gene set of interest (for example, the differentailly expressed genes) in cell type1/cell type2 for the calculation of the nichenetr based ligand/receptor activity scores.  
#' @param lr_network: Ligand-receptor network matrix with two columns "from" and "to", consisting of the gene symbols of the candidate ligands and receptors, respectively, in pairs.
#' @param ligand_target_matrix_ct1_to_ct2: Matrix of ligand to target regulatory potential scores calculated by nichenetr's algorithm, from ligands expressed in cell type1 to targets in cell type2.  Ligand gene symbols in columns and target gene symbols in rows.
#' @param ligand_target_matrix_ct2_to_ct1: Matrix of ligand to target regulatory potential scores calculated by nichenetr's algorithm, from ligands expressed in cell type2 to targets in cell type1.  Ligand gene symbols in columns and target gene symbols in rows.
#' @param receptor_target_matrix_ct1_to_ct2: Matrix of receptor to target regulatory potential scores calculated by nichenetr's algorithm, from receptors expressed in cell type2 (that are cognate receptors of ligands expressed in cell type1) to targets in cell type2.  Receptor gene symbols in columns and target gene symbols in rows.
#' @param receptor_target_matrix_ct2_to_ct1: Matrix of receptor to target regulatory potential scores calculated by nichenetr's algorithm, from receptors expressed in cell type1 (that are cognate receptors of ligands expressed in cell type2) to targets in cell type1.  Receptor gene symbols in columns and target gene symbols in rows.
#' @param discrete_error_rate: Value of the error_rate variable of the nichenetr package function "make_discrete_ligand_target_matrix" -- FDR for cutoff_method "fdrtool" and "distribution"; number between 0 and 1 indicating which top fraction of target genes should be returned for cutoff_method "quantile".
#' @param discrete_cutoff_method: Value of the cutoff_method variable of the nichenetr package function "make_discrete_ligand_target_matrix" -- Method to determine which genes can be considered as a target and which genes not, based on the regulatory probability scores. Possible options: "distribution", "fdrtool" and "quantile".
#' @param discrete_fdr_method: Value of the fdr_method variable of the nichenetr package function "make_discrete_ligand_target_matrix" -- Only relevant when cutoff_method is "fdrtool". Possible options: "global" and "local"
#' @return
#'  Basics: A list with the following elements:
#'  Basics$ave_expr_ct1: A matrix of average gene expression values (LogNormalized) in each condition in cell type1.  Gene symbols in rows, conditions in columns
#'  Basics$ave_expr_ct2: A matrix of average gene expression values (LogNormalized) in each condition in cell type2.  Gene symbols in rows, conditions in columns
#'  Basics$pct_expr_ct1: A matrix of gene detection rates in each condition in cell type1.  Gene symbols in rows, conditions in columns
#'  Basics$pct_expr_ct2: A matrix of gene detection rates in each condition in cell type2.  Gene symbols in rows, conditions in columns
#'  Basics$thresh_expr_ct1: A binary matrix of 0s and 1s indicating if the detection rate of each gene passed the specified cutoff value min_pct (1: Yes; 0: No) in each condition in cell type1.  Gene symbols in rows, conditions in columns
#'  Basics$thresh_expr_ct2: A binary matrix of 0s and 1s indicating if the detection rate of each gene passed the specified cutoff value min_pct (1: Yes; 0: No) in each condition in cell type2.  Gene symbols in rows, conditions in columns
#'  Basics$genes_thresh_expr_ct1: A vector of gene symbols that are "expressed" (passed detection rate cutoff min_pct) in at least one condition in cell type1
#'  Basics$genes_thresh_expr_ct2: A vector of gene symbols that are "expressed" (passed detection rate cutoff min_pct) in at least one condition in cell type2
#'  Basics$lr_expr_ct1_to_ct2: A list with two elements:
#'  Basics$lr_expr_ct1_to_ct2$eachcondition: A list of ligand (expressed in cell type1)-receptor (expressed in cell type2) pairs "expressed" in each condition
#'  Basics$lr_expr_ct1_to_ct2$bind: Ligand (expressed in cell type1)-receptor (expressed in cell type2) pairs "expressed" in at least one condition
#'  Basics$lr_expr_ct2_to_ct1: A list with two elements:
#'  Basics$lr_expr_ct2_to_ct1$eachcondition: A list of ligand (expressed in cell type2)-receptor (expressed in cell type1) pairs "expressed" in each condition
#'  Basics$lr_expr_ct2_to_ct1$bind: Ligand (expressed in cell type2)-receptor (expressed in cell type1) pairs "expressed" in at least one condition
#'  Basics$ligand_activities_matrix_ct1_to_ct2: Matrix of regulatory potential scores of candidate ligands from cell type1 (calcuated by nichenetr's algorithm)
#'  Basics$ligand_activities_matrix_ct2_to_ct1: Matrix of regulatory potential scores of candidate ligands from cell type2 (calcuated by nichenetr's algorithm)
#'  Basics$receptor_activities_matrix_ct1_to_ct2: Matrix of regulatory potential scores of candidate receptors in cell type2 (that are cognate receptors of candidate ligands from cell type1) (calcuated by nichenetr's algorithm)
#'  Basics$receptor_activities_matrix_ct2_to_ct1: Matrix of regulatory potential scores of candidate receptors in cell type1 (that are cognate receptors of candidate ligands from cell type2) (calcuated by nichenetr's algorithm)
#'  Basics$myLRL: List of identified LRLoops:
#'  Basics$myLRL$R1->L2_R2->L1: Matrix of LRLoops with columns "L1", "R1", "L2" and "R2" where L1-R1 are candidate ligand-receptor pairs from cell type1 to cell type2, L2-R2 are candidate ligand-receptor pairs from cell type2 to cell type1; L2 is a top target of R1, L1 is a top target of R2
#'  Basics$myLRL$L1->L2_L2->L1: Matrix of LRLoops with columns "L1", "R1", "L2" and "R2" where L1-R1 are candidate ligand-receptor pairs from cell type1 to cell type2, L2-R2 are candidate ligand-receptor pairs from cell type2 to cell type1; L2 is a top target of L1, L1 is a top target of L2
#'
#'
#' @import nichenetr tidyverse Seurat
#' @export


PrepareBasics <- function(ct1obj, ct2obj, min_pct, geneset_ct1, geneset_ct2,
                          lr_network, ligand_target_matrix_ct1_to_ct2, receptor_target_matrix_ct1_to_ct2, ligand_target_matrix_ct2_to_ct1, receptor_target_matrix_ct2_to_ct1,discrete_error_rate, discrete_cutoff_method, discrete_fdr_method) {
  ######
  library(Seurat)
  library(Matrix)
  
  ######
  conditions = as.vector(unique(ct1obj@meta.data[,'Condition']))
  
  ### Expression data (LogNormalized)
  data_ct1 = ct1obj@assays$RNA@data
  data_ct2 = ct2obj@assays$RNA@data
  
  ### Average gene expression values in each condition 
  ave_expr_ct1 = log1p(AverageExpression(ct1obj, group.by = 'Condition')$RNA)
  ave_expr_ct2 = log1p(AverageExpression(ct2obj, group.by = 'Condition')$RNA)
  if (length(conditions) > 1) {
    ave_expr_ct1 = ave_expr_ct1[,conditions]
    ave_expr_ct2 = ave_expr_ct2[,conditions]
  } else {
    colnames(ave_expr_ct1) = conditions
    colnames(ave_expr_ct2) = conditions
  }
  
  ### Gene detection rates
  print('1')
  pct_expr_ct1 = get_pct_expr(ct1obj, conditions)
  pct_expr_ct2 = get_pct_expr(ct2obj, conditions)
  
  ### Identify the "expressed" genes with a detection rate greater than a particular threshold
  print('2')
  thresh_expr_ct1 = (pct_expr_ct1 > min_pct)*1
  thresh_expr_ct2 = (pct_expr_ct2 > min_pct)*1
  genes_thresh_expr_ct1 = rownames(thresh_expr_ct1)[rowSums(thresh_expr_ct1)>0]
  genes_thresh_expr_ct2 = rownames(thresh_expr_ct2)[rowSums(thresh_expr_ct2)>0]
  
  ### Identify "expressed" ligand-receptor pairs
  lr_expr_ct1_to_ct2 = get_expr_lr(lr_network = lr_network,
                                   thresh_expr_from = thresh_expr_ct1, thresh_expr_to = thresh_expr_ct2, conditions = conditions)
  lr_expr_ct2_to_ct1 = get_expr_lr(lr_network = lr_network, 
                                   thresh_expr_from = thresh_expr_ct2, thresh_expr_to = thresh_expr_ct1, conditions = conditions)
  
  ### Print number of ligand-recepor pairs from ct1 to ct2 and from ct2 to ct1
  numlr = matrix(0, nrow = (1+length(conditions)), ncol = 2)
  rownames(numlr) = c(conditions, 'bind')
  colnames(numlr) = c('num_lr_expr_ct1_to_ct2', 'num_lr_expr_ct2_to_ct1')
  for (i in 1:length(conditions)) {
    numlr[i,'num_lr_expr_ct1_to_ct2'] = nrow(lr_expr_ct1_to_ct2[['eachcondition']][[conditions[i]]])
    numlr[i,'num_lr_expr_ct2_to_ct1'] = nrow(lr_expr_ct2_to_ct1[['eachcondition']][[conditions[i]]])
  }
  numlr['bind', 'num_lr_expr_ct1_to_ct2'] = nrow(lr_expr_ct1_to_ct2[['bind']])
  numlr['bind', 'num_lr_expr_ct2_to_ct1'] = nrow(lr_expr_ct2_to_ct1[['bind']])
  print("number of <expressed> ligand-receptor pairs:")
  print(numlr)
  
  if (sum(numlr==0)>0) 
    stop("Please relax the criteria defining <expressed> genes, or remove the condition(s) with no <expressed> ligand-receptor pairs in either direction (ct1 to ct2, or, ct2 to ct1).")
  
  ### Ligand and receptor regulatory potential analysis (nichenetr scores)
  background_expr_genes_ct2 = genes_thresh_expr_ct2 %>% .[. %in% rownames(ligand_target_matrix_ct1_to_ct2)]
  ligand_activities_ct1_to_ct2 = predict_ligand_activities(geneset = geneset_ct2, 
                                                           background_expressed_genes = background_expr_genes_ct2, 
                                                           ligand_target_matrix = ligand_target_matrix_ct1_to_ct2, 
                                                           potential_ligands = intersect(unique(lr_expr_ct1_to_ct2$bind[,'from']),
                                                                                         colnames(ligand_target_matrix_ct1_to_ct2)))
  ligand_activities_matrix_ct1_to_ct2 = as.matrix(ligand_activities_ct1_to_ct2 %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson))))
  background_expr_genes_ct2 = genes_thresh_expr_ct2 %>% .[. %in% rownames(receptor_target_matrix_ct1_to_ct2)]
  receptor_activities_ct1_to_ct2 = predict_ligand_activities(geneset = geneset_ct2,
                                                             background_expressed_genes = background_expr_genes_ct2, 
                                                             ligand_target_matrix = receptor_target_matrix_ct1_to_ct2, 
                                                             potential_ligands = intersect(unique(lr_expr_ct1_to_ct2$bind[,'to']), 
                                                                                           colnames(receptor_target_matrix_ct1_to_ct2)))
  receptor_activities_matrix_ct1_to_ct2 = as.matrix(receptor_activities_ct1_to_ct2 %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson))))
  colnames(receptor_activities_matrix_ct1_to_ct2)[1] = 'test_receptor'
  
  background_expr_genes_ct1 = genes_thresh_expr_ct1 %>% .[. %in% rownames(ligand_target_matrix_ct2_to_ct1)]
  ligand_activities_ct2_to_ct1 = predict_ligand_activities(geneset = geneset_ct1, 
                                                           background_expressed_genes = background_expr_genes_ct1, 
                                                           ligand_target_matrix = ligand_target_matrix_ct2_to_ct1, 
                                                           potential_ligands = intersect(unique(lr_expr_ct2_to_ct1$bind[,'from']),
                                                                                         colnames(ligand_target_matrix_ct2_to_ct1)))
  ligand_activities_matrix_ct2_to_ct1 = as.matrix(ligand_activities_ct2_to_ct1 %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson))))
  background_expr_genes_ct1 = genes_thresh_expr_ct1 %>% .[. %in% rownames(receptor_target_matrix_ct2_to_ct1)]
  receptor_activities_ct2_to_ct1 = predict_ligand_activities(geneset = geneset_ct1,
                                                             background_expressed_genes = background_expr_genes_ct1, 
                                                             ligand_target_matrix = receptor_target_matrix_ct2_to_ct1, 
                                                             potential_ligands = intersect(unique(lr_expr_ct2_to_ct1$bind[,'to']), 
                                                                                           colnames(receptor_target_matrix_ct2_to_ct1)))
  receptor_activities_matrix_ct2_to_ct1 = as.matrix(receptor_activities_ct2_to_ct1 %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson))))
  colnames(receptor_activities_matrix_ct2_to_ct1)[1] = 'test_receptor'
  
  ### Identify the set of targets of each ligand/receptor
  ligand_target_matrix_binary_ct2_to_ct1 = make_discrete_ligand_target_matrix(ligand_target_matrix_ct2_to_ct1,
                                                                              error_rate = discrete_error_rate, cutoff_method = discrete_cutoff_method,
                                                                              fdr_method = discrete_fdr_method,ligands_position = "cols")*1
  ligand_target_matrix_binary_ct1_to_ct2 = make_discrete_ligand_target_matrix(ligand_target_matrix_ct1_to_ct2, 
                                                                              error_rate = discrete_error_rate, cutoff_method = discrete_cutoff_method, 
                                                                              fdr_method = discrete_fdr_method,ligands_position = "cols")*1
  receptor_target_matrix_binary_ct2_to_ct1 = make_discrete_ligand_target_matrix(receptor_target_matrix_ct2_to_ct1, 
                                                                                error_rate = discrete_error_rate, cutoff_method = discrete_cutoff_method, 
                                                                                fdr_method = discrete_fdr_method,ligands_position = "cols")*1
  receptor_target_matrix_binary_ct1_to_ct2 = make_discrete_ligand_target_matrix(receptor_target_matrix_ct1_to_ct2, 
                                                                                error_rate = discrete_error_rate, cutoff_method = discrete_cutoff_method, 
                                                                                fdr_method = discrete_fdr_method,ligands_position = "cols")*1
  ### Get LRloop network
  myLRL = get_LRL(lr_expr_ct2_to_ct1, lr_expr_ct1_to_ct2, 
                  ligand_target_matrix_binary_ct2_to_ct1, ligand_target_matrix_binary_ct1_to_ct2, 
                  receptor_target_matrix_binary_ct2_to_ct1, receptor_target_matrix_binary_ct1_to_ct2)
  
  
  
  
  ### Collect basics 
  Basics = list(ave_expr_ct1, ave_expr_ct2, pct_expr_ct1, pct_expr_ct2, thresh_expr_ct1, thresh_expr_ct2, genes_thresh_expr_ct1, genes_thresh_expr_ct2,
                lr_expr_ct1_to_ct2, lr_expr_ct2_to_ct1,
                ligand_activities_matrix_ct1_to_ct2, receptor_activities_matrix_ct1_to_ct2, 
                ligand_activities_matrix_ct2_to_ct1, receptor_activities_matrix_ct2_to_ct1,
                myLRL)
  names(Basics) = c("ave_expr_ct1", "ave_expr_ct2", "pct_expr_ct1", "pct_expr_ct2", 
                    "thresh_expr_ct1", "thresh_expr_ct2", "genes_thresh_expr_ct1", "genes_thresh_expr_ct2",
                    "lr_expr_ct1_to_ct2", "lr_expr_ct2_to_ct1",
                    "ligand_activities_matrix_ct1_to_ct2", "receptor_activities_matrix_ct1_to_ct2", 
                    "ligand_activities_matrix_ct2_to_ct1", "receptor_activities_matrix_ct2_to_ct1",
                    "myLRL")
  
  return(Basics)
  
}
