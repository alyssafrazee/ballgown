\name{assmb2annot}
\alias{assmb2annot}
\title{match assembled transcript(s) to closest annotated transcript/gene}
\usage{
assmb2annot(gtf, chr, assembled)
}
\arguments{
  \item{gtf}{path to a GTF file containing locations of
  annotated transcripts}

  \item{chr}{which chromosome is your assembled transcript
  of interest on? Should be formatted the same way as
  chromosomes in \code{gtf}.}

  \item{assembled}{GRanges object, with ranges representing
  the assembled transcript's exons.}
}
\value{
a GRangesList containing the annotated transcript(s) that
most overlap \code{assembled}.
}
\description{
match assembled transcript(s) to closest annotated
transcript/gene
}
\author{
Alyssa Frazee
}

