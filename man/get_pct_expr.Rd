% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_pct_expr.R
\name{get_pct_expr}
\alias{get_pct_expr}
\title{Calculate the gene detection rate in each condition}
\usage{
get_pct_expr(seuratobj, conditions)
}
\arguments{
\item{seuratobj}{A Seurat object}

\item{conditions}{A vector of conditions of interest}
}
\value{
A matrix of the detection rate of each gene in each condition with genes in rows and conditions in columns
}
\description{
Calculate the gene detection rate in each condition
}
