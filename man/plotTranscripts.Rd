\name{plotTranscripts}
\alias{plotTranscripts}
\title{visualize structure of assembled transcripts}
\usage{
plotTranscripts(gene, gown, samples = NULL, colorby = "transcript",
  meas = "FPKM", legend = TRUE, labelTranscripts = FALSE, main = NULL,
  colorBorders = FALSE, log = FALSE, logbase = 2)
}
\arguments{
  \item{gene}{name of gene whose transcripts will be
  plotted.  When using Cufflinks output, usually of the
  form \code{"XLOC_######"}}

  \item{gown}{ballgown object containing experimental and
  phenotype data}

  \item{samples}{vector of sample(s) to plot. Can be
  \code{'none'} if only one plot (showing transcript
  structure in gray) is desired. Use
  \code{sampleNames(gown)} to see sample names for
  \code{gown}. Defaults to \code{sampleNames(gown)[1]}.}

  \item{colorby}{one of \code{"transcript"}, \code{"exon"},
  or \code{"none"}, indicating which feature's abundances
  should dictate plot coloring.  If \code{"none"}, all
  transcripts are drawn in gray.}

  \item{meas}{which expression measurement to color
  features by, if any. Must match an available measurement
  for whatever feature you're plotting.}

  \item{legend}{if \code{TRUE} (as it is by default), a
  color legend is drawn on top of the plot indicating
  scales for feature abundances.}

  \item{labelTranscripts}{if \code{TRUE}, transcript ids
  are labeled on the left side of the plot. Default
  \code{FALSE}.}

  \item{main}{optional string giving the desired plot
  title.}

  \item{colorBorders}{if \code{TRUE}, exon borders are also
  drawn in color (instead of black, as they are by
  default). Useful for visualizing abundances for skinny
  exons and/or smaller plots, as often happens when
  \code{length(samples)} is large.}

  \item{log}{if \code{TRUE}, color transcripts on the log
  scale. Default \code{FALSE}. To account for expression
  values of 0, we add 1 to all expression values before
  taking the log.}

  \item{logbase}{log base to use if \code{log = TRUE}.
  default 2.}
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

