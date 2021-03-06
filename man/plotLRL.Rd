% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotLRL.R
\name{plotLRL}
\alias{plotLRL}
\title{Plot the LRloop network 
(A_ and B_ in the LR node lables: Distinguish ligand-from-ct1_receptor-in-ct2 pairs from ligand-from-ct2_receptor-in-ct1 pairs)}
\usage{
plotLRL(
  LRloop_info,
  WhichL1R1L2R2Score,
  nodecolors,
  labelcolor,
  nodesize,
  edgecutoff,
  labelsize,
  edgecolor,
  pheatmapcolor,
  pheatmap_dist
)
}
\arguments{
\item{LRloop_info}{The list of the LRloop_network info calculated by function "LRL_info_collection"}

\item{WhichL1R1L2R2Score}{Run the function "print_LRL_Score_options" first to check available options}

\item{nodecolors}{Colors of the A_L1_R1 and B_L2_R2 type nodes in a vector of length two}

\item{labelcolor}{Label text color}

\item{nodesize}{Node size}

\item{edgecutoff}{Only plot edges with L1R1L2R2 LRloopScore greater than edgecutoff}

\item{labelsize}{Node label size}

\item{edgecolor}{Edge color}

\item{pheatmapcolor}{Value of pheatmap variable 'color'}

\item{pheatmap_dist}{Value of pheatmap variables 'clustering_distance_rows' and 'clustering_distance_cols'. e.g. 'correlation', 'euclidean'}
}
\description{
Plot the LRloop network 
(A_ and B_ in the LR node lables: Distinguish ligand-from-ct1_receptor-in-ct2 pairs from ligand-from-ct2_receptor-in-ct1 pairs)
}
