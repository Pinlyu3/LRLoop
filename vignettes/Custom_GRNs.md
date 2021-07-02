# Perform LRLoop analysis using custom GRNs
 
  - A. Construct ligand_target_matrix as in https://github.com/saeyslab/nichenetr/blob/master/vignettes/model_construction.md
  - B. Similar to the construction of ligand_target_matrix, with some functions in nichenetr modified slightly, construct receptor_target_matrix
  - As an example, in our work, we combine the ligand-receptor networks from NicheNet and connectomeDB2020
  The signaling network and gene regulatory network are taken directly from NicheNet
  Inputs in this example: 
  NicheNet_lr_network.rds
  NATMI_lrc2p.csv
  NicheNet_signaling_network.rds
  NicheNet_gr_network.rds
  Outputs in this example (needed for further analysis):
  lr_network
  ligand_target_matrix
  receptor_target_matrix

``` r
library(LRLoop)
```

## Step 1: prepair the ligand-receptor, signaling and gene regulatory networks

### ligand-receptor network:

``` r
lr_network_nichenet = readRDS("Networks/NicheNet_lr_network.rds") # The ligand-receptor network collected in NicheNet
lr_network_NATMI = read.csv("Networks/NATMI_lrc2p.csv") # The ligand-receptor network collected in connectomeDB2020 for NATMI (Network Analysis Toolkit for the Multicellular Interactions)
lr_network_NATMI = cbind(lr_network_NATMI, rep('NATMI', nrow(lr_network_NATMI)), rep('NATMI', nrow(lr_network_NATMI)))
colnames(lr_network_NATMI) = c('from', 'to', 'source', 'database')
lr_network_NATMI = tibble(lr_network_NATMI)
lr_network_large = rbind(lr_network_nichenet, lr_network_NATMI) # Combine the ligand-receptor networks of interest
# Remark: as another example, we use the pre-filtered literature supported lr_network with the ligands and recepros filtered by the annotations in NATMI and CelltalkDB:
# lr_network = read.table("Networks/LR_network_202103.txt", header = TRUE)
# lr_network_large = tibble(lr_network[,c(1:4)])
``` 

### Any non-NicheNet-collected networks should be added to the model:

``` r
new_network_weights_df = tibble(source = 'NATMI', weight = 1)
source_weights_df = source_weights_df %>% bind_rows(new_network_weights_df)
new_annotation_data_source = tibble(source = 'NATMI', database = 'NATMI', type_db = 'NATMI', 
                                    type_interaction = 'NATMI', network = 'ligand_receptor')
annotation_data_sources = rbind(annotation_data_sources, new_annotation_data_source)

### signaling network 
sig_network = readRDS("Networks/NicheNet_signaling_network.rds")

### gene regulatory network 
gr_network = readRDS("Networks/NicheNet_gr_network.rds")

### if user provided context-specific gene regulatory network (named "mygr_network") with source and database "mygr" is preferred, then after "mygr_network" is loaded and processed (optional):
 gr_network = mygr_network
 new_network_weights_df = tibble(source = 'mygr', weight = 1)
 source_weights_df = source_weights_df %>% bind_rows(new_network_weights_df)
 new_annotation_data_source = tibble(source = 'mygr', database = 'mygr', type_db = 'mygr', 
                                     type_interaction = 'mygr', network = 'gene_regulatory')
 annotation_data_sources = rbind(annotation_data_sources, new_annotation_data_source)
# Remark: If mygr_network is prefered, also remove the original gr_network's data sources from source_weights_df in the next step

### Only keep selected data sources (optional): 
# In this example, we remove the ppi predicted ligand-receptor pairs from the ligand-receptor network
data_sources_to_remove = annotation_data_sources %>% filter(database %in% c("ppi_prediction","ppi_prediction_go")) %>% pull(source)
data_sources_to_keep = annotation_data_sources$source %>% setdiff(data_sources_to_remove) 
source_weights_df = source_weights_df %>% filter(source %in% data_sources_to_keep)
lr_network = lr_network_large %>% filter(source %in% data_sources_to_keep) 
``` 

## Step 2-A: Calculate the ligand_target_matrix by NicheNet

### Aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network

``` r
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, 
source_weights_df = source_weights_df)
``` 

### Downweight the importance of signaling and gene regulatory hubs 

``` r
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, 
lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)
```


### Calculate target gene regulatory potential scores for specified ligands

``` r
ligands = unique(as.list(as.data.frame(lr_network)[,'from']))
ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, 
                                                      algorithm = "PPR", 
                                                      damping_factor = hyperparameter_list$damping_factor, 
                                                      ltf_cutoff = hyperparameter_list$ltf_cutoff)
```



## Step 2-B: Calculate the receptor_target_matrix by NicheNet's algorithm


### Aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
``` r
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, 
source_weights_df = source_weights_df, n_output_networks = 3)
```
### Downweight the importance of signaling and gene regulatory hubs
``` r
weighted_networks = apply_hub_corrections_noligand(weighted_networks = weighted_networks, 
sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)
```

### Calculate target gene regulatory potential scores for specified receptors
``` r
signaling_network = weighted_networks$sig
regulatory_network = weighted_networks$gr
allgenes =  unique(c(signaling_network$from, signaling_network$to, regulatory_network$from, regulatory_network$to))
receptors = unique(as.list(intersect(as.data.frame(lr_network)[,'to'], allgenes)))
receptor_target_matrix = construct_receptor_target_matrix(weighted_networks = weighted_networks, receptors = receptors, 
                                                          damping_factor = hyperparameter_list$damping_factor, 
                                                          rtf_cutoff = hyperparameter_list$ltf_cutoff)
``` 

## Convert to mouse genes, if needed

``` r
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

colnames(receptor_target_matrix) = receptor_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(receptor_target_matrix) = receptor_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
receptor_target_matrix = receptor_target_matrix %>% .[!is.na(rownames(receptor_target_matrix)), !is.na(colnames(receptor_target_matrix))]
``` 
