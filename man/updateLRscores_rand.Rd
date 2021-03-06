% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updateLRscores_rand.R
\name{updateLRscores_rand}
\alias{updateLRscores_rand}
\title{Get modified LRscores with random LRloop networks}
\usage{
updateLRscores_rand(
  n,
  lr_ct1_to_ct2,
  lr_ct2_to_ct1,
  LRloops,
  LRS,
  LRC,
  from_to,
  LRCform,
  lambda,
  mu,
  k
)
}
\arguments{
\item{n}{Number of random LRloop networks to be constructed.}

\item{lr_ct1_to_ct2}{A matrix with columns "from" and "to" of ligand-receptor pairs from ct1 to ct2.}

\item{lr_ct2_to_ct1}{A matrix with columns "from" and "to" of ligand-receptor pairs from ct2 to ct1.}

\item{LRloops}{The original LRloop network of all possible LRloops btween ct1 and ct2 of the LRpairs in lr_network.}

\item{LRS}{Same inputs as those used in updateLRscores.}

\item{LRC}{Same inputs as those used in updateLRscores.}

\item{from_to}{Same inputs as those used in updateLRscores.}

\item{LRCform}{Same inputs as those used in updateLRscores.}

\item{lambda}{Same inputs as those used in updateLRscores.}

\item{mu}{Same inputs as those used in updateLRscores.}

\item{k}{Same inputs as those used in updateLRscores.}
}
\value{
LRscoresUpdated_nrand: Updated LRscore with n random LRloop networks.
}
\description{
Get modified LRscores with random LRloop networks
}
