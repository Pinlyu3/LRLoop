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

## LRloop analysis
### Identify the conditions of interest 
``` r
conditions = unique(ct1obj@meta.data[,'Condition'])

#> conditions
#[1] "mmP60"      "mmNMDA03"   "mmNMDA06"   "mmNMDA12"   "mmNMDA24"   "mmNMDA36"   "mmNMDA48"   "mmNMDA48FI" "mmNMDA72"  
```
### Differential expression analysis for all condition pairs of interest (if length(conditions) >= 2)
 Remark: Users can use other methods to find differentially expressed genes of interest, just note that the results of this step should be two lists "DEGinfo_ct1" and "DEGinfo_ct2" in which:
 "DEGinfo_ct1$DEG" and "DEGinfo_ct2$DEG" are lists where each element of them stores the DEA results for one pair of conditions of interest, which is a matrix with genes in rows and at least two colomns "ave_log2FC" and "p_val_adj";  
 "DEGinfo_ct1$DEgenes" and "DEGinfo_ct2$DEgenes" are vectors of differentially expressed gene symbols (the union of the DEGs resulted from all pairs of conditions of interest)
 In this example, for each cell type (microglia and MG), we perform DEA between each NMDA time point and mmP60 (control) and collect the union of all the DEGs
 ``` r
DEGinfo_ct1 = get_DEG(seuratobj = ct1obj, idents_1 = conditions[2:9], idents_2 = conditions[1], 
                      only_pos = FALSE, min_pct = 0.1, logfc_threshold = 0.25, p_val_adj_threshold = 0.05, test_use = "wilcox")
DEGinfo_ct2 = get_DEG(seuratobj = ct2obj, idents_1 = conditions[2:9], idents_2 = conditions[1], 
                      only_pos = FALSE, min_pct = 0.1, logfc_threshold = 0.25, p_val_adj_threshold = 0.05, test_use = "wilcox")
                    
#> names(DEGinfo_ct1)
#[1] "DEG"     "DEgenes"
#> names(DEGinfo_ct1$DEG)
#[1] "mmNMDA03_vs_mmP60"   "mmNMDA06_vs_mmP60"   "mmNMDA12_vs_mmP60"   "mmNMDA24_vs_mmP60"   "mmNMDA36_vs_mmP60"   "mmNMDA48_vs_mmP60"   "mmNMDA48FI_vs_mmP60" "mmNMDA72_vs_mmP60"  
#> head(DEGinfo_ct1$DEG$mmNMDA03_vs_mmP60)
#             p_val avg_log2FC pct.1 pct.2    p_val_adj
#Fth1  4.077072e-21   3.664240 1.000 0.793 1.138849e-16
#Srgn  3.063184e-20   3.410460 0.969 0.103 8.556391e-16
#Ftl1  2.086013e-18   2.617280 1.000 0.931 5.826861e-14
#Hmox1 3.108638e-18   5.565917 0.862 0.052 8.683357e-14
#Ccl3  2.758505e-17   4.892835 0.908 0.224 7.705331e-13
#Id2   2.901389e-17   4.377565 0.815 0.034 8.104450e-13
#> head(DEGinfo_ct1$DEgenes)
#[1] "Fth1"  "Srgn"  "Ftl1"  "Hmox1" "Ccl3"  "Id2"  
```


