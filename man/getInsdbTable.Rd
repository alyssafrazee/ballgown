\name{getInsdbTable}
\alias{getInsdbTable}
\title{get feature-by-sample expression table from cuffdiff2/cummeRbund output}
\usage{
getInsdbTable(dbFile, feature = "isoform", meas = "fpkm")
}
\arguments{
  \item{dbFile}{path to the database file from cummeRbund
  (available as an insilicodb download)}

  \item{feature}{one of \code{"isoform"}, \code{"gene"},
  \code{"CDS"}, or \code{"TSS"} -- what genomic feature do
  you want expression measurements for? Default
  \code{"isoform"}}

  \item{meas}{one of \code{"fpkm"}, \code{"raw_frags"},
  \code{"internal_scaled_frags"},
  \code{"external_scaled_frags"} -- which measurement from
  the Cuffdiff output table do you want? Default
  \code{"fpkm"}.}
}
\value{
matrix with \code{feature}s in rows, samples (replicates)
in columns, and \code{meas}urements in cells
}
\description{
get feature-by-sample expression table from
cuffdiff2/cummeRbund output
}
\details{
on a test dataset with 37K transcripts and 24 reps, this
took about 30 seconds
}
\author{
Alyssa Frazee
}

