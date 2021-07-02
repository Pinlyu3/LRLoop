# Perform LRLoop analysis starting form a Seurat object

## load  packages:
``` r
library(LRLoop)
```
## load data:
### Required inputs:
-  Seurat object (data LogNormalized) of celltype1 (ct1) and celltype2 (ct2) "ct1obj" and "ct2obj":
- a) Their corresponding expression data "ct1obj@assays$RNA@data" and "ct2obj@assays$RNA@data" should have the same rows (genes in rows); 
- b) Their metadata "ct1obj@meta.data" and "ct2obj@meta.data" both have the column "Condition" with the same set of condtions of interest.
-  lr_network with columns "from" (ligands in this column) and "to" (receptors in this column).
-  ligand_target_matrix_ct1_to_ct2 & ligand_target_matrix_ct2_to_ct1: resulted from Construct_ligandreceptor_target_matrix.R. Users can also choose to use a default pre-calculated ligand_target_matrix.
-  receptor_target_matrix_ct1_to_ct2 & receptor_target_matrix_ct2_to_ct1: resulted from Construct_ligandreceptor_target_matrix.R. Users can also choose to use a default pre-calculated receptor_target_matrix.


Load the Seurat objects "ct1obj" and "ct2obj" for celltype1 and celltype2 of interest
``` r
load("ExampleData/ct1obj.RData")
load("ExampleData/ct2obj.RData")

#> ct1obj
#An object of class Seurat 
#27933 features across 631 samples within 1 assay 
#Active assay: RNA (27933 features, 0 variable features)
#> ct2obj
#An object of class Seurat 
#27933 features across 5942 samples within 1 assay 
#Active assay: RNA (27933 features, 0 variable features)
#> unique(ct1obj@meta.data[,'Condition'])
#[1] "mmP60"      "mmNMDA03"   "mmNMDA06"   "mmNMDA12"   "mmNMDA24"   "mmNMDA36"   "mmNMDA48"   "mmNMDA48FI" "mmNMDA72"  
#> unique(ct2obj@meta.data[,'Condition'])
#[1] "mmP60"      "mmNMDA03"   "mmNMDA06"   "mmNMDA12"   "mmNMDA24"   "mmNMDA36"   "mmNMDA48"   "mmNMDA48FI" "mmNMDA72"
``` 
Load ligand-receptor network and ligand/receptor-target potential matrices
Remark: In this example, these ligand/receptor_target_matrix_* were calculated based on the lr_network, sig_network and gr_networks filtered by the genes detected in microglia or MG, respectively
``` r
load("ExampleData/lr_network.RData")

#> head(lr_network)
## A tibble: 6 x 4
#  from   to    source         database
#  <chr>  <chr> <chr>          <chr>   
#1 Cxcl1  Cxcr2 kegg_cytokines kegg    
#2 Cxcl2  Cxcr2 kegg_cytokines kegg    
#3 Cxcl5  Cxcr2 kegg_cytokines kegg    
#4 Ppbp   Cxcr2 kegg_cytokines kegg    
#5 Cxcl9  Cxcr3 kegg_cytokines kegg    
#6 Cxcl10 Cxcr3 kegg_cytokines kegg    

lr_network = unique(as.matrix(lr_network)[,c('from','to')]) 
lr_network = lr_network[lr_network[,'from'] %in% rownames(ct1obj@assays$RNA@data) & lr_network[,'to'] %in% rownames(ct1obj@assays$RNA@data),]

#> head(lr_network)
#     from     to     
#[1,] "Cxcl1"  "Cxcr2"
#[2,] "Cxcl2"  "Cxcr2"
#[3,] "Cxcl5"  "Cxcr2"
#[4,] "Ppbp"   "Cxcr2"
#[5,] "Cxcl9"  "Cxcr3"
#[6,] "Cxcl10" "Cxcr3"

load("ExampleData/ligand_target_matrix_ct1_to_ct2.RData")
load("ExampleData/ligand_target_matrix_ct2_to_ct1.RData")
load("ExampleData/receptor_target_matrix_ct1_to_ct2.RData")
load("ExampleData/receptor_target_matrix_ct2_to_ct1.RData")

#> ligand_target_matrix_ct1_to_ct2[1:5,1:5]
#        Cxcl1        Cxcl2        Cxcl5         Ppbp       Cxcl12
#A1bg    0.0008333858 0.0008145524 0.0008655993 0.0006291958 0.0011008142
#A1cf    0.0008510756 0.0005766328 0.0005491118 0.0004963418 0.0007971723
#A2m     0.0014920135 0.0010681160 0.0010995387 0.0011364595 0.0020007205
#A3galt2 0.0005118905 0.0002699805 0.0002650488 0.0002536862 0.0003670101
#A4galt  0.0012909658 0.0008870256 0.0009401015 0.0008000234 0.0014363265

# Or just use a set of our pre-calculated default networks:
#load("Networks/DefaultNetworks/mm/NicheNet_and_NATMI_bonafide_filtered/lr_network.RData")
#load("Networks/DefaultNetworks/mm/NicheNet_and_NATMI_bonafide_filtered/ligand_target_matrix.RData")
#load("Networks/DefaultNetworks/mm/NicheNet_and_NATMI_bonafide_filtered/receptor_target_matrix.RData")
# Then run:
#lr_network = unique(as.matrix(lr_network)[,c('from','to')]) 
#lr_network = lr_network[lr_network[,'from'] %in% rownames(ct1obj@assays$RNA@data) & lr_network[,'to'] %in% rownames(ct1obj@assays$RNA@data),]
#ligand_target_matrix_ct1_to_ct2 = ligand_target_matrix
#ligand_target_matrix_ct2_to_ct1 = ligand_target_matrix
#receptor_target_matrix_ct1_to_ct2 = receptor_target_matrix
#receptor_target_matrix_ct2_to_ct1 = receptor_target_matrix
``` 
