% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/construct_tf_target_matrix_noligand.R
\name{construct_tf_target_matrix_noligand}
\alias{construct_tf_target_matrix_noligand}
\title{Convert gene regulatory network into a matrix.}
\usage{
construct_tf_target_matrix_noligand(
  weighted_networks,
  tfs_as_cols = FALSE,
  standalone_output = FALSE
)
}
\arguments{
\item{weighted_networks}{the Seurat object of interest}

\item{tfs_as_cols}{vectors of idents}

\item{standalone_output}{vectors of idents}
}
\value{
A matrix of tf-target regulatory weights
}
\description{
Convert gene regulatory network into a matrix.
}
