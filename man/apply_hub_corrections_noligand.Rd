% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply_hub_corrections_noligand.R
\name{apply_hub_corrections_noligand}
\alias{apply_hub_corrections_noligand}
\title{Downweight the hub nodes in the signaling and gene regulatory networks}
\usage{
apply_hub_corrections_noligand(weighted_networks, sig_hub, gr_hub)
}
\arguments{
\item{weighted_networks}{A list of three weighted networks lr, sig and gr (ligand-receptor network, signaling network and gene regulatory network), in data frame/tibble with columns "from", "to" and "weight"}

\item{sig_hub}{A number between 0 and 1. 0: no correction for hubiness; 1: maximal correction for hubiness.}

\item{gr_hub}{A number between 0 and 1. 0: no correction for hubiness; 1: maximal correction for hubiness.}
}
\value{
A list of two hubiness-corrected weighted networks sig and gr (signaling network and gene regulatory network), in data frame/tibble with columns "from", "to" and "weight".
}
\description{
Downweight the hub nodes in the signaling and gene regulatory networks
}
