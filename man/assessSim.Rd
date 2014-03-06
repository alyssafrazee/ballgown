\name{assessSim}
\alias{assessSim}
\title{analyze results of a two-group simulation experiment}
\usage{
assessSim(bg, bgresults, annotation, chr, trulyDEids, cuffdiffFile,
  qcut = 0.05, UCSC = TRUE, ret = FALSE, nClosest = 1)
}
\arguments{
  \item{bg}{ballgown object from the simulated data}

  \item{bgresults}{data frame resulting from a call to
  \code{stattest(bg,...)}. Should be a transcript-level
  test (i.e., \code{feature="transcript"} in
  \code{stattest}).}

  \item{annotation}{either a path to an annotation gtf
  file, or a data.frame of an annotation file that was
  already read in and subset to exons on the right
  chromosome.}

  \item{chr}{}

  \item{trulyDEids}{}

  \item{cuffdiffFile}{}

  \item{qcut}{either a number between 0 and 1 to be used as
  the q-value significance cutoff, or a vector like
  \code{seq(0,1,by=0.01)} (i.e., ranging from 0 to 1 in
  even increments).}

  \item{UCSC}{}

  \item{nClosest}{for true positives, make a positive call
  if any of the \code{nClosest} closest (most overlapping)
  transcripts are called DE by a method. Default 1.}
}
\value{
creates ROC curve and returns sensitivities/specificities
used for said ROC curve, comparing cuffdiff to ballgown (as
used to create \code{bg}).
}
\description{
[documentation in progress]
}
\details{
\code{trulyDEids} should be the transcripts that were set
to be differentially expressed, identified the SAME WAY AS
THEY ARE in \code{annotation}. This is super important!!!

Also: if \code{qcut} is a vector, \code{assessSim} returns
sensitivities/specificities and creates an ROC plot. If
\code{qcut} is a single number, sensitivity and specificity
for both methods are returned, using \code{qcut} as a
q-value significance cutoff.
}
\author{
Alyssa Frazee
}

