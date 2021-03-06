% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRL_info_collection.R
\name{LRL_info_collection}
\alias{LRL_info_collection}
\title{Collect data on LRloop network relevant ligands and receptors}
\usage{
LRL_info_collection(
  valuse_ct1,
  valuse_ct2,
  scalar,
  ScoreConditionIdx,
  LRscore_method,
  LoopScore_method,
  LRlfcscore_method,
  LooplfcScore_method,
  LRloop_network,
  LRL_eachcondition,
  ave_expr_ct1,
  ave_expr_ct2,
  pct_expr_ct1,
  pct_expr_ct2,
  DEGinfo_ct1,
  DEGinfo_ct2,
  genes_cluster_ct1,
  genes_cluster_ct2,
  L1R1_cluster,
  L2R2_cluster,
  ligand_activities_matrix_ct1_to_ct2,
  receptor_activities_matrix_ct1_to_ct2,
  ligand_activities_matrix_ct2_to_ct1,
  receptor_activities_matrix_ct2_to_ct1
)
}
\arguments{
\item{valuse_ct1}{Matrix of the values to be used in the calculation of LRscores and L1R1L2R2 LoopScores with genes in rows and conditions in columns corresponding to ct1/ct2.}

\item{valuse_ct2}{Matrix of the values to be used in the calculation of LRscores and L1R1L2R2 LoopScores with genes in rows and conditions in columns corresponding to ct1/ct2.}

\item{scalar}{A number, plays the role of a scalar in the calculation of LRscores when LRscore_method is set to 'scsigr' or individual_scale'.}

\item{ScoreConditionIdx}{Column indices of valuse_ct1 (same as that of valuse_ct2) to be used in the calculation of the maximum or mean LRscores and LoopScores across conditions}

\item{LRscore_method}{Method used in the calculation of LRscores, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr'}

\item{LoopScore_method}{Method used in the calculation of L1R1L2R2 LoopScores, available options are "ave_geo" and "ave_alg"}

\item{LRlfcscore_method}{Method used in the calculation of logFC-based LRscores, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr'}

\item{LooplfcScore_method}{Method used in the calculation of logFC-based L1R1L2R2 LoopScores, available options are "ave_geo" and "ave_alg"}

\item{LRloop_network}{An LRloop network matrix with columns "L1", "R1", "L2" and "R2".}

\item{LRL_eachcondition}{List of LRloop networks in each condition.}

\item{ave_expr_ct1}{Matrix of average expression value of each gene in each condition in ct1/ct2 cells with genes in rows and conditions in columns.}

\item{ave_expr_ct2}{Matrix of average expression value of each gene in each condition in ct1/ct2 cells with genes in rows and conditions in columns.}

\item{pct_expr_ct1}{Matrix of fraction of ct1/ct2 cells that expressed each gene in each condition with genes in rows and conditions in columns.}

\item{pct_expr_ct2}{Matrix of fraction of ct1/ct2 cells that expressed each gene in each condition with genes in rows and conditions in columns.}

\item{DEGinfo_ct1}{List of two elements.  The first element DEGinfo_ct1/ct2$DEG is a list where each element is a matrix with genes in rows and at least two colomns "ave_log2FC" and "p_val_adj";  The second element DEGinfo_ct1/ct2$DEgenes is a vector that contains the gene symbols of DEGs in ct1/ct2.}

\item{DEGinfo_ct2}{List of two elements.  The first element DEGinfo_ct1/ct2$DEG is a list where each element is a matrix with genes in rows and at least two colomns "ave_log2FC" and "p_val_adj";  The second element DEGinfo_ct1/ct2$DEgenes is a vector that contains the gene symbols of DEGs in ct1/ct2.}

\item{genes_cluster_ct1}{User defined ct1/ct2 gene clustering vector (integers 0,1,2,...) with the name of each element the corresponding gene symbol}

\item{genes_cluster_ct2}{User defined ct1/ct2 gene clustering vector (integers 0,1,2,...) with the name of each element the corresponding gene symbol}

\item{L1R1_cluster}{User defined ct1 to ct2/ct2 to ct1 ligand-receptor pair clustering vector (integers 0,1,2,...) with the name of each element the corresponding ligand-receptor gene symbols (in the form of L-R)}

\item{L2R2_cluster}{User defined ct1 to ct2/ct2 to ct1 ligand-receptor pair clustering vector (integers 0,1,2,...) with the name of each element the corresponding ligand-receptor gene symbols (in the form of L-R)}

\item{ligand_activities_matrix_ct1_to_ct2}{Ligand/receptor_activities_matrix_ct1_to_ct2 calculated by nichenetr's algorithm}

\item{receptor_activities_matrix_ct1_to_ct2}{Ligand/receptor_activities_matrix_ct1_to_ct2 calculated by nichenetr's algorithm}

\item{ligand_activities_matrix_ct2_to_ct1}{Ligand/receptor_activities_matrix_ct2_to_ct1 calculated by nichenetr's algorithm}

\item{receptor_activities_matrix_ct2_to_ct1}{Ligand/receptor_activities_matrix_ct2_to_ct1 calculated by nichenetr's algorithm}
}
\value{
LRL_info_collection: a list of matrices:
 LRL_info_collection$'LRloop expression': The LRloop network, together with the average expression value of each gene and the fraction of cells expressing each gene in the corresponding cell type and each condition
 LRL_info_collection$'L1R1 expression': Identdified ligand-receptor pairs from ct1 to ct2, together with the average expression value of each gene and the fraction of cells expressing each gene in the corresponding cell type and each condition
 LRL_info_collection$'L2R2 expression': Identified ligand-receptor pairs from ct2 to ct1, together with the average expression value of each gene and the fraction of cells expressing each gene in the corresponding cell type and each condition
 LRL_info_collection$'LRloopDEG': Differential expresssion information of the genes in LRloop_network 
 LRL_info_collection$'L1R1DEG': Differential expression information of L1R1 pairs
 LRL_info_collection$'L2R2DEG': Differential expression information of L2R2 pairs
 LRL_info_collection$'L1R1L2R2 clusters': User defined clusters of L1, R1, L2, R2, L1-R1 pair and L2-R2 pair
 LRL_info_collection$'LRloop expr_score': The L1R1L2R2 LoopScores calculated based on valuse_ct1 and valuse_ct2
 LRL_info_collection$'LRloop logFC_score': The L1R1L2R2 LoopScores calculated based on the absolute values of the logFC of L1, R1, L2 and R2
 LRL_info_collection$'LRloop nichenet_score': The nichenet algorithm based scores of L1, R1, L2 and L2 in each L1R1L2R2 loop
 LRL_info_collection$'L1R1 expr_score': The LRscores of L1R1 pairs calculated based on valuse_ct1 and valuse_ct2
 LRL_info_collection$'L1R1 logFC_score': The LRscores of L1R1 pairs calculated based on the absolute values of the logFC of L1 and R1
 LRL_info_collection$'L1R1 nichenet_score': The nichenet algorithm based scores of L1 and R1 in each L1R1 pair
 LRL_info_collection$'L2R2 expr_score': The LRscores of the L2R2 pairs calculated based on valuse_ct2 and valuse_ct1
 LRL_info_collection$'L2R2 logFC_score': The LRscores of the L2R2 pairs calculated based on the absolute values of the logFC of L2 and R2
 LRL_info_collection$'L2R2 nichenet_score': The nichenet algorithm based scores of L2 and R2 in each L2R2 pair
 LRL_info_collection$'L1R1 expr_score-LRLoop_filtered': The LRscores of L1R1 pairs calculated based on valuse_ct1 and valuse_ct2 and filtered by LRloops in each condition, that is, in each condition, an LRscore will be set to 0 if the L1R1 pair does not form any LRloop in that conditon. 
 LRL_info_collection$'L2R2 expr_score-LRLoop_filtered': The LRscores of L2R2 pairs calculated based on valuse_ct2 and valuse_ct1 and filtered by LRloops in each condition, that is, in each condition, an LRscore will be set to 0 if the L2R2 pair does not form any LRloop in that conditon.
}
\description{
Collect data on LRloop network relevant ligands and receptors
}
