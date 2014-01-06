\docType{methods}
\name{structure}
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
\alias{structure}
\alias{subset}
\alias{texpr}
\usage{
structure(x)

data(x)

indexes(x)

dirs(x)

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

  \item{x}{ballgown object}

  \item{x}{ballgown object return list containing elements
  \code{e2t}, \code{i2t}, \code{t2g}, \code{bamfiles}, and
  \code{pData}, where \code{e2t} and \code{i2t} are data
  frames linking exons and introns (respectively) to
  transcripts, \code{t2g} is a data frame linking
  transcripts to genes, and \code{bamfiles} and
  \code{pData} are described at the \code{link{ballgown}}
  constructor help page.}

  \item{x}{ballgown object}

  \item{x}{ballgown object}

  \item{...}{either \code{'cov'}, \code{'FPKM'}, depending
  on which type of expression measurement you want to
  extract (provide no argument to extract both).}

  \item{...}{\code{'rcount'}, \code{'ucount'},
  \code{'mrcount'}, \code{'cov'}, \code{'cov_sd'},
  \code{'mcov'}, or \code{'mcov_sd'}, depending on which
  type of expression measurement you want to extract
  (provide no argument to extract all).}

  \item{...}{\code{'rcount'}, \code{'ucount'}, or
  \code{'mrcount'}, depending on which type of expression
  measurement you want to extract (provide no argument to
  extract all).}
}
\value{
list containing elements \code{intron} (GRanges),
\code{exon} (GRanges), and \code{trans} (GRangesList),
denoting genomic positions of exons, introns, and
transcripts (represented as sets of exons).

list containing elements \code{intron}, \code{exon}, and
\code{trans} (all data frames) -- feature-by-sample
expression tables.

paths to the on-disk directories holding the data (created
with \code{tablemaker}) used to create \code{x}

the date \code{x} was created
}

