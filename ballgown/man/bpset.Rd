\name{bpset}
\alias{bpset}
\title{get vector of nucleotides comprising a transcript}
\usage{
  bpset(tra)
}
\arguments{
  \item{tra}{a GRanges object representing a transcript
  (see the transcript "structure" element of a ballgown
  object)}
}
\value{
  vector of nucleotide positions included in the transcript
}
\description{
  get vector of nucleotides comprising a transcript
}
\seealso{
  \code{\link{transcriptOverlaps}} which is this function's
  wrapper
}

