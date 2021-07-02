# Perform LRLoop analysis starting form a Seurat object

## load  packages:
``` r
library(LRLoop)
```
## load data
### Required inputs:
-  Seurat object (data LogNormalized) of celltype1 (ct1) and celltype2 (ct2) "ct1obj" and "ct2obj":
- a) Their corresponding expression data "ct1obj@assays$RNA@data" and "ct2obj@assays$RNA@data" should have the same rows (genes in rows); 
- b) Their metadata "ct1obj@meta.data" and "ct2obj@meta.data" both have the column "Condition" with the same set of condtions of interest.
-  lr_network with columns "from" (ligands in this column) and "to" (receptors in this column).
-  ligand_target_matrix_ct1_to_ct2 & ligand_target_matrix_ct2_to_ct1: resulted from Construct_ligandreceptor_target_matrix.R. Users can also choose to use a default pre-calculated ligand_target_matrix.
-  receptor_target_matrix_ct1_to_ct2 & receptor_target_matrix_ct2_to_ct1: resulted from Construct_ligandreceptor_target_matrix.R. Users can also choose to use a default pre-calculated receptor_target_matrix.
