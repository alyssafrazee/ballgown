\name{gffReadGR}
\alias{gffReadGR}
\title{read in gtf file as GRanges object}
\usage{
gffReadGR(gtf)
}
\arguments{
  \item{gtf}{name of GTF/GFF file on disk}
}
\value{
object of class \code{GRanges} representing the genomic
features in \code{gtf}
}
\description{
(very) light wrapper for rtracklayer::import
}
\author{
Alyssa Frazee
}
\seealso{
\code{\link{gffRead}} for reading in a GTF file as a data
frame rather than GRanges.
}

