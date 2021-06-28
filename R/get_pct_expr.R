# Calculate the gene detection rate in each condition
# Inputs: 
# seuratobj: the Seurat object of interest
# data: normalized expressed data
# conditions: conditions/time points of interest
# Output:
# pct_expr: a matrix with genes in rows and conditions in columns and each entry the percentage of cells expressing each gene in each condition

get_pct_expr <- function(seuratobj, data, conditions) {
  pct_expr = matrix(0, nrow = nrow(data), ncol = length(conditions))
  Idents(object = seuratobj) = 'Condition'
  for (i in 1:length(conditions)) {
    subdata = data[,WhichCells(seuratobj, idents = conditions[i])]
    pct_expr[,i] = rowSums(subdata >0)/ncol(subdata)
  }
  rownames(pct_expr) = rownames(data)
  colnames(pct_expr) = conditions
  
  return(pct_expr)
}




