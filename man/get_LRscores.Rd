% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_LRscores.R
\name{get_LRscores}
\alias{get_LRscores}
\title{Calculate the LRscore of each "expressed" ligand-receptor pair in each condition.  Only run after running "PrepareBasics".}
\usage{
get_LRscores(
  lr_expr,
  conditions,
  value.use_from,
  value.use_to,
  scalar,
  LRscore_method,
  thresh_expr_from,
  thresh_expr_to,
  LRL_eachcondition,
  LRL_filter = "none",
  thresh_expr_cut = FALSE
)
}
\arguments{
\item{lr_expr}{A ligand-receptor network matrix with two columns "from" and "to".}

\item{conditions}{A vector of conditions of interest}

\item{value.use_from}{A numeric matrix with genes in rows and conditions in columns, selected as the ligand values to be used in the calculation of LRscores}

\item{value.use_to}{A numeric matrix with genes in rows and conditions in columns, selected as the receptor values to be used in the calculation of LRscores}

\item{scalar}{A number, plays the role of a scalar in the calculation of LRscores when the LRscore_method set to 'scsigr' or individual_scale'}

\item{LRscore_method}{The method of calculating the LRscores, available options are 'mean', 'individual_scale', 'individual_scale_exp', 'product', 'bias_receptor' and 'scsigr'}

\item{thresh_expr_from}{A numeric matrix of 0s and 1s with genes in rows and conditions in columns (with the same column names as value.use_from) that define whether each gene is expressed in each condition (0-NO, 1-YES) in the sender cells.}

\item{thresh_expr_to}{A numeric matrix of 0s and 1s with genes in rows and conditions in columns (with the same column names as value.use_to) that define whether each gene is expressed in each condition (0-NO, 1-YES) in the receiver cells.}

\item{LRL_eachcondition}{List of LRloop networks in each condition.}

\item{LRL_filter}{"L1R1", "L2R2" or "none".  "L1R1" ("L2R2"): The LR pairs are from ct1 (ct2) to ct2 (ct1) and apply the LRloop filter in each condition, that is, in each condition, an LRscore will be set to 0 if the LR pair does not form any LRloop in that conditon. "none": No LRloop filter applied.}

\item{thresh_expr_cut}{TRUE or FALSE.  If TRUE, the LRscore of ligand-receptor pairs for which either the ligand or the recepor is not expressed as defined in thresh_expr_from and thresh_expr_to, respectively, is set to 0.  Default is FALSE.}
}
\value{
A matrix with columns "ligand", "receptor" and the condition names that keeps record of the LRscore of each ligand-receptor pair in each condition
}
\description{
Calculate the LRscore of each "expressed" ligand-receptor pair in each condition.  Only run after running "PrepareBasics".
}
