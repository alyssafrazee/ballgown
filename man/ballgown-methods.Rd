\docType{methods}
\name{ballgown-methods}
\alias{ballgown-methods}
\alias{data}
\alias{data<-}
\alias{dirs}
\alias{eexpr}
\alias{gexpr}
\alias{iexpr}
\alias{indexes}
\alias{indexes<-}
\alias{mergedDate}
\alias{pData}
\alias{pData<-}
\alias{sampleNames}
\alias{structure}
\alias{subset}
\alias{texpr}
\title{methods for objects of class \code{ballgown}}
\usage{
structure(x)

data(x)

indexes(x)

dirs(x)

sampleNames(x)

mergedDate(x)

indexes(x) <- value

data(x) <- value

subset(x, ...)

pData(x)

pData(x) <- value

texpr(x, ...)

eexpr(x, ...)

iexpr(x, ...)

gexpr(x)
}
\arguments{
  \item{x}{ballgown object}

  \item{value}{the updated value for a ballgown object
  component}

  \item{...}{for \code{subset}: logical expression
  involving any column of the \code{texpr} dictating what
  subset of \code{x} is desired.  For \code{*expr} methods:
  one of \code{'cov'}, \code{'FPKM'}, \code{'rcount'},
  \code{'ucount'}, \code{'mrcount'}, \code{'cov_sd'},
  \code{'mcov'}, or \code{'mcov_sd'}, depending on which
  type of expression measurement is desired.  Leave
  \code{...} blank to select all expression measurements.}

  \item{meas}{expression measurement to extract (for use
  with \code{*expr} methods).  Defaults to FPKM for
  \code{texpr} and \code{rcount} for \code{eexpr} and
  \code{iexpr}. Specifying \code{"all"} will return all
  expression measurements as well as extra feature-level
  information.}
}
\value{
for \code{structure}: list containing elements
\code{intron} (GRanges), \code{exon} (GRanges), and
\code{trans} (GRangesList), denoting genomic positions of
exons, introns, and transcripts (represented as sets of
exons).

for \code{data}: list containing elements \code{intron},
\code{exon}, and \code{trans} (all data frames) --
feature-by-sample expression tables.

for \code{indexes}: list containing elements \code{e2t},
\code{i2t}, \code{t2g}, \code{bamfiles}, and \code{pData},
where \code{e2t} and \code{i2t} are data frames linking
exons and introns (respectively) to transcripts, \code{t2g}
is a data frame linking transcripts to genes, and
\code{bamfiles} and \code{pData} are described at the
\code{link{ballgown}} constructor help page.

for \code{dirs}: paths to the on-disk directories holding
the data (created with \code{tablemaker}) used to create
\code{x}

for \code{sampleNames}: names of samples included in the
dataset.  Matches folder names of \code{dirs(x)}.

for \code{mergedDate}: the date \code{x} was created

for \code{*expr} methods: a feature-by-sample table with
the specified expression measurement(s) in the cells.
}
\description{
The following methods are available for the S4 class
\code{ballgown}.
}

