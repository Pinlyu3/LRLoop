# Identify the "expressed" ligand-receptor pairs
# Inputs: 
# lr_network: the ligand-receptor network
# thresh_expr_from: the binary matrix indicating whether each gene is "expressed" in each condition in the sender cells in this step
# thresh_expr_to: the binary matrix indicating whether each gene is "expressed" in each condition in the receiver cells in this step
# conditions: conditions/time points of interest
# Output:
# lr_expr: a list:
# lr_expr$eachcondition: a list, ligand-receptor pairs "expressed" in each condition
# lr_expr$bind: ligand-receptor pairs that "expressed" in at least one condition

get_expr_lr <- function(lr_network, thresh_expr_from, thresh_expr_to, conditions) {
  
  lr_expr_list = list()
  for (i in 1:length(conditions)) {
    lr_expr_list[[i]] = lr_network[thresh_expr_from[lr_network[,'from'], conditions[i]]==1 & thresh_expr_to[lr_network[,'to'], conditions[i]]==1,]
  }
  names(lr_expr_list) = conditions
  
  lr_expr_bind = vector()
  for (i in 1:length(lr_expr_list)) {
    lr_expr_bind = rbind(lr_expr_bind, lr_expr_list[[i]])
  }
  lr_expr_bind = unique(lr_expr_bind)
  
  lr_expr = list(lr_expr_list, lr_expr_bind)
  names(lr_expr) = c('eachcondition', 'bind')
  
  return(lr_expr)
}


