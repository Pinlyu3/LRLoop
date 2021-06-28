# Calculate the LRscore of each "expressed" ligand-receptor pair in each condition
# Inputs: 
# lr_expr: the "expressed" ligand-receptor pairs
# conditions: the conditions of interest
# value.use_from: ave_expr_sender or pct_expr_sender
# value.use_to: ave_expr_receiver or pct_expr_receiver
# scalar: a number 
# LRscore_method: the method of calculating the LRscores, one of 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr' 
# Output:
# LRscorematrix: a matrix with columns "ligand", "receptor" and the condition names that keeps record of the LRscore of each ligand-receptor pair in each condition 

get_LRscores <- function(lr_expr, conditions, value.use_from, value.use_to, scalar, LRscore_method) {
  
  LRscorematrix = matrix(0, nrow = nrow(lr_expr), ncol = (2+length(conditions)))
  colnames(LRscorematrix) = c('from', 'to', conditions)
  for (i in 1:nrow(LRscorematrix)) {
    ligand = lr_expr[i,'from']
    receptor = lr_expr[i,'to']
    lrscores = vector()
    for (j in 1:length(conditions)) {
      lrscores[j] = LRscore(value.use_from[ligand,conditions[j]], value.use_to[receptor,conditions[j]], scalar, LRscore_method)
    }
    LRscorematrix[i,] = c(ligand, receptor, lrscores)
  }
  
  return(LRscorematrix)
}



