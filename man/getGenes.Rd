\name{getGenes}
\alias{getGenes}
\title{label assembled transcripts with gene names}
\usage{
getGenes(gtf, assembled, UCSC = TRUE)
}
\arguments{
  \item{gtf}{path to a GTF file containing locations of
  annotated transcripts}

  \item{assembled}{GRangesList object, with each set of
  ranges representing exons of an assembled transcript.}

  \item{UCSC}{set to \code{TRUE} if you're using a UCSC gtf
  file. (Requires some extra text-processing).}
}
\value{
a character vector of the same length as \code{assembled},
providing the name(s) of the gene(s) that overlaps each
transcript in \code{assembled}.
}
\description{
label assembled transcripts with gene names
}
\details{
chromosome labels in \code{gtf} and \code{assembled} should
match. (i.e., you should provide the path to a gtf
corrsponding to the same annotation you used when
constructing \code{assembled})
}
\author{
Alyssa Frazee
}

