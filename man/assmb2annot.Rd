\name{assmb2annot}
\alias{assmb2annot}
\title{match assembled transcript(s) to closest annotated transcript/gene}
\usage{
assmb2annot(gtf, assembled)
}
\arguments{
  \item{gtf}{path to a GTF file containing locations of
  annotated transcripts}

  \item{chr}{which chromosome is your assembled transcript
  of interest on? Should be formatted the same way as
  chromosomes in \code{gtf}.}

  \item{assembled}{GRangesList object, with each set of
  ranges representing exons of an assembled transcript.}
}
\value{
a GRangesList containing the corresponding annotated
transcript(s) that most overlap \code{assembled}.
}
\description{
match assembled transcript(s) to closest annotated
transcript/gene
}
\details{
the \code{elementMetadata} slot of each \code{GRanges}
object in the returned \code{GRangesList} contains the
annotated \code{gene_id} and \code{transcript_id}. The
\code{names} component of this \code{GRangesList} gives the
assembled transcript ids.

Also be careful not to confuse this with
\code{annot2assmb}, which finds the closest
\emph{assembled} transcript to each \emph{annotated}
transcript. That function is more useful in simulations;
this one is more useful for getting biological information
out of an assembly.
}
\author{
Alyssa Frazee
}

