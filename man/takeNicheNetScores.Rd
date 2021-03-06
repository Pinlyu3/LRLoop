% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/takeNicheNetScores.R
\name{takeNicheNetScores}
\alias{takeNicheNetScores}
\title{Put ligand/receptor_activity scores for expressed ligand-receptor pairs from NicheNet in a specific format that can be further used in the calculation of LRLoop corrected scores.}
\usage{
takeNicheNetScores(
  LRpairs,
  from_to,
  Basics,
  nichenetscore = "pearson",
  method = "Ligand"
)
}
\arguments{
\item{LRpairs}{A matrix of LR pairs of interest with columns "from" and "to", in which each row is an LR pair with ligand in the column "from" and receptor in the column "to".}

\item{from_to}{Specify the LRpairs are ligands from celltype1(ct1) to receptors in celltype2(ct2) or ligands from celltype2(ct2) to receptors in celltype1(ct1).}

\item{Basics}{Result from the function PrepareBasics}

\item{nichenetscore}{NicheNet score options, "auroc", "aupr" or "pearson".  Default is "pearson".}

\item{method}{Options are "Ligand", "Receptor", "LRmean_arithmetic" or "LRmean_geometric".  Default is "Ligand". If method = "Ligand", for each LRpair, only use the ligand activity score. If method = "Receptor", for each LRpair, only use the receptor activity score. If method = "LRmean_arithmetic", for each LRpair, use the arithmetic average of the ligand activity score and the receptor activity score. If method = "LRmean_geometric", for each LRpair, use the geometric average of the ligand activity score and the receptor activity score.}
}
\value{
A matrix with three columns "from", "to" and "score" with the ligand, receptor and the score value in each column, respectively.
}
\description{
Put ligand/receptor_activity scores for expressed ligand-receptor pairs from NicheNet in a specific format that can be further used in the calculation of LRLoop corrected scores.
}
