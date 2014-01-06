\name{plotMeans}
\alias{plotMeans}
\title{visualize transcript abundance by group}
\usage{
plotMeans(gene, gown, overall = TRUE, groupvar = "none",
  groupname = "none", dattype = c("cov", "FPKM", "rcount", "ucount",
  "mrcount", "mcov"), colorby = c("transcript", "exon"), legend = TRUE,
  labelTranscripts = FALSE)
}
\arguments{
  \item{gene}{name of gene whose transcripts will be
  plotted.  When using Cufflinks output, usually of the
  form \code{"XLOC_######"}}

  \item{gown}{ballgown object containing experimental and
  phenotype data}

  \item{overall}{if \code{TRUE}, color features by the
  overall (experiment-wide) mean rather than a
  group-specific mean}

  \item{groupvar}{string representing the name of the
  variable denoting which sample belongs to which group.
  Can be \code{"none"} (if you want the study-wide mean),
  or must correspond to the name of a column of
  \code{pData(gown)}.  Usually a categorical variable.}

  \item{groupname}{string representing which group's
  expression means you want to plot.  Can be \code{"none"}
  (if you want the study-wide mean), \code{"all"} (if you
  want a multipanel plot of each group's mean expression),
  or any of the levels of \code{groupvar}.}

  \item{colorby}{one of \code{"transcript"} or
  \code{"exon"}, indicating which feature's abundances
  should dictate plot coloring.}

  \item{legend}{if \code{TRUE} (as it is by default), a
  color legend is drawn on top of the plot indicating the
  scale for feature abundances.}

  \item{labelTranscripts}{if \code{TRUE}, transcript ids
  are labeled on the left side of the plot. Default
  \code{FALSE}.}
}
\value{
produces a plot of the transcript structure for the
specified gene in the current graphics device, colored by
study-wide or group-specific mean expression level.
}
\description{
visualize transcript abundance by group
}
\author{
Alyssa Frazee
}
\seealso{
\code{\link{plotTranscripts}}
}

