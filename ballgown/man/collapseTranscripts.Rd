\name{collapseTranscripts}
\alias{collapseTranscripts}
\title{cluster a gene's transcripts and calculate cluster-level expression}
\usage{
collapseTranscripts(gown, dattype = c("cov", "FPKM"), method = c("hclust",
  "kmeans"))
}
\arguments{
  \item{gown}{ballgown object containing experimental data}

  \item{dattype}{which transcript-level expression
  measurement to use (\code{'cov'}, average per-base
  coverage, or \code{'FPKM'})}

  \item{method}{which clustering method to use
  (\code{'hclust'}, hierarchical clustering, or
  \code{'kmeans'}, k-means clustering)}
}
\value{
data frame with one row per transcript cluster and one
column per sample, where entries are summed expression
measurements for all the transcripts in the appropriate
cluster
}
\description{
cluster a gene's transcripts and calculate cluster-level
expression
}
\details{
Transcript clustering methods are in development, so use
this function with caution: for example, it's not clear
that the appropriate cluster-level expression measurement
is the sum.

Also, this function runs clustering and collapsing on the
entire ballgown object and could be very slow, so you may
want to check out the \code{\link{subset}} method for
ballgown objects and run this function on small chunks of
genes.
}
\author{
Alyssa Frazee
}
\seealso{
\code{\link{hclust}}, \code{\link{\kmeans}},
\code{\link{clusterTranscripts}} for gene-level transcript
clustering, \code{\link{plotLatentTranscripts}} for
visualizing transcript clusters
}

