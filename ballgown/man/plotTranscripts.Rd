\name{plotTranscripts}
\alias{plotTranscripts}
\title{visualize structure of assembled transcripts}
\usage{
plotTranscripts(gene, samp = NULL, gown, legend = TRUE,
  labelTranscripts = FALSE, colorby = c("transcript", "exon", "none"),
  main = NULL)
}
\arguments{
  \item{gene}{name of gene whose transcripts will be
  plotted.  When using Cufflinks output, usually of the
  form \code{"XLOC_######"}}

  \item{samp}{name of sample whose abundances dictate
  transcript colors.  Can be \code{NULL} if no colors are
  desired. Usually prefixed with the transcript expression
  measurement followed by a dot, e.g.,
  \code{"cov.sample1"}}

  \item{gown}{ballgown object containing experimental and
  phenotype data}

  \item{legend}{if \code{TRUE} (as it is by default), a
  color legend is drawn on top of the plot indicating
  scales for feature abundances.}

  \item{labelTranscripts}{if \code{TRUE}, transcript ids
  are labeled on the left side of the plot. Default
  \code{FALSE}.}

  \item{colorby}{one of \code{"transcript"}, \code{"exon"},
  or \code{"none"}, indicating which feature's abundances
  should dictate plot coloring.  If \code{"none"}, all
  transcripts are drawn in gray.}

  \item{main}{optional string giving the desired plot
  title.}
}
\value{
produces a plot of the transcript structure for the
specified gene in the current graphics device.
}
\description{
visualize structure of assembled transcripts
}
\author{
Alyssa Frazee
}
\seealso{
\code{\link{plotMeans}}
}

