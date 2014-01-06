\name{checkAssembledTx}
\alias{checkAssembledTx}
\title{plot annotated and assembled transcripts together}
\usage{
checkAssembledTx(tx.assembled, tx.annotated, ind = 1,
  main = "Assembled and Annotated Transcripts")
}
\arguments{
  \item{tx.assembled}{a GRangesList object where the
  GRanges objects in the list represent sets of exons
  comprising assembled transcripts}

  \item{tx.annotated}{a GRangesList object where the
  GRanges objects in the list represent sets of exons
  comprising annotated transcripts}

  \item{ind}{integer; index of \code{tx.annotated}
  specifying which annotated transcript to plot.  All
  transcripts (assembled and annotated) overlapping
  \code{tx.annotated[[ind]]} will be plotted.  Default 1.}

  \item{main}{optional character string giving the title
  for the resulting plot.  Default: "Assembled and
  Annotated Transcripts"}
}
\value{
No return value, but a plot is produced with annotated
transcripts on the bottom panel (shaded in gray) and
assembled transcripts on the top panel (shaded with
diagonal lines).
}
\description{
plot annotated and assembled transcripts together
}
\author{
Alyssa Frazee
}

