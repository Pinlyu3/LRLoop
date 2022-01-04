#' Calculate some basic statistics on LRloops
#' @name LRL_stats
#' @param LRL_list List of identified LRloops calculated by the function "get_LRL"
#' @import nichenetr tidyverse Seurat stats
#' @export


LRL_stats <- function(LRL_list) {
  slots = names(LRL_list)[2:length(names(LRL_list))]
  p_A_a_list_stats = matrix(0, nrow = length(slots), ncol = 15)
  colnames(p_A_a_list_stats) = c('reg_type', 
                                 'num_[(L1_R1)-(L2_R2)]_pairs', 'num_L1R1_pairs', 'num_L2R2_pairs', 'num_L1', 'num_R1', 'num_L2', 'num_R2',
                                 'mean_num(L2_R2)_each(L1_R1)', 'median_num(L2_R2)_each(L1_R1)', 'num_[(L1_R1)-(L2_R2)]_pairs_unordered',
                                 'num_[(L1)-(L2)]_pairs', 'mean_num(L2)_each(L1)', 'median_num(L2)_each(L1)', 'num_[(L1)-(L2)]_pairs_unordered')
  rownames(p_A_a_list_stats) = slots
  for (k in 1:length(slots)) {
    A = LRL_list[[slots[k]]]
    p_A_a_list_stats[k,'reg_type'] = slots[k]
    p_A_a_list_stats[k,'num_[(L1_R1)-(L2_R2)]_pairs'] = nrow(A)
    p_A_a_list_stats[k,'num_L1'] = length(unique(A[,'L1']))
    p_A_a_list_stats[k,'num_R1'] = length(unique(A[,'R1']))
    p_A_a_list_stats[k,'num_L2'] = length(unique(A[,'L2']))
    p_A_a_list_stats[k,'num_R2'] = length(unique(A[,'R2']))
    
    A = cbind(sprintf("%s__%s", A[,'L2'], A[,'R2']), sprintf("%s__%s", A[,'L1'], A[,'R1']))
    colnames(A) = c('L2__R2', 'L1__R1')
    p_A_a_list_stats[k,'num_L1R1_pairs'] = length(unique(A[,'L1__R1']))
    p_A_a_list_stats[k,'num_L2R2_pairs'] = length(unique(A[,'L2__R2']))
    
    LRs = unique(A[,'L1__R1'])
    numaA = vector()
    for(i in 1:length(LRs)) {
      numaA[i] = sum(A[,'L1__R1']==LRs[i])
    }
    p_A_a_list_stats[k,'mean_num(L2_R2)_each(L1_R1)'] = mean(numaA)
    p_A_a_list_stats[k,'median_num(L2_R2)_each(L1_R1)'] = median(numaA)
    
    A = A[!duplicated(t(apply(A, 1, sort))),]
    p_A_a_list_stats[k,'num_[(L1_R1)-(L2_R2)]_pairs_unordered'] = nrow(A)
    
    if(slots[k] %in% c('L1->L2 & L2->L1')) {
      A = unique(LRL_list[[slots[k]]][,c('L2', 'L1')])
      p_A_a_list_stats[k,'num_[(L1)-(L2)]_pairs'] = nrow(A)
      
      Ls = unique(A[,'L1'])
      numa = vector()
      for(i in 1:length(Ls)) {
        numa[i] = sum(A[,'L1']==Ls[i])
      }
      p_A_a_list_stats[k,'mean_num(L2)_each(L1)'] = mean(numa)
      p_A_a_list_stats[k,'median_num(L2)_each(L1)'] = median(numa)
      
      A = A[!duplicated(t(apply(A, 1, sort))),]
      p_A_a_list_stats[k,'num_[(L1)-(L2)]_pairs_unordered'] = nrow(A)
    }
  }
  
  return(p_A_a_list_stats)
}

