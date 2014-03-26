\name{get_read_numbers}
\alias{get_read_numbers}
\title{get read numbers from zero-inflated negative binomial simulation model}
\usage{
get_read_numbers(fc, nbparams, modfit, libsizes, seed = NULL)
}
\arguments{
  \item{fc}{matrix of fold changes. Should have number of
  rows equal to the number of transcripts in the
  simulation, and number of columns equal to the number of
  samples in the simulation.}

  \item{nbparams}{zero-inflated negative binomial
  parameters (estimated with \code{get_params})}

  \item{modfit}{model defining the mean/dispersion
  relationship for transcript expression. Usually estimated
  with \code{get_params}.}

  \item{libsizes}{how many reads should be generated for
  each sample? Should have length equal to number of
  columns of \code{fc}, i.e., \code{libsizes} is a vector
  with length equal to the number of samples in the
  simulation.}

  \item{seed}{optional seed (for reproducibility)}
}
\description{
get read numbers from zero-inflated negative binomial
simulation model
}
\author{
Jeff Leek
}

