\name{gffReadGR}
\alias{gffReadGR}
\title{read in gtf file as GRanges object}
\usage{
gffReadGR(gtf, splitByTranscript = FALSE, identifier = "transcript_id",
  sep = ";")
}
\arguments{
  \item{gtf}{name of GTF/GFF file on disk}
}
\value{
if \code{splitByTranscript} is \code{FALSE}, an object of
class \code{GRanges} representing the genomic features in
\code{gtf}. If \code{splitByTranscript} is TRUE, an object
of class \code{GRangesList}, where each element is a
\code{GRanges} object corresponding to an annotated
transcript (designated in \code{names}).
}
\description{
(very) light wrapper for rtracklayer::import
}
\author{
Alyssa Frazee
}
\seealso{
\code{\link{gffRead}} for reading in a GTF file as a data
frame rather than a \code{GRanges}/\code{GRangesList}
object.
}

