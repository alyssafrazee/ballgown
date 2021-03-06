% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGenes.R
\name{getGenes}
\alias{getGenes}
\title{label assembled transcripts with gene names}
\usage{
getGenes(gtf, assembled, UCSC = TRUE, attribute = "gene_id")
}
\arguments{
\item{gtf}{path to a GTF file containing locations of annotated transcripts}

\item{assembled}{GRangesList object, with each set of ranges representing 
exons of an assembled transcript.}

\item{UCSC}{set to \code{TRUE} if you're using a UCSC gtf file. (Requires
some extra text processing).}

\item{attribute}{set to attribute name in \code{gtf} that gives desired gene
identifiers. Default \code{"gene_id"}; another commone one is
\code{"gene_name"} (for the gene symbol).}
}
\value{
an \code{IRanges CharacterList} of the same length as
  \code{assembled}, providing the name(s) of the gene(s) that overlaps each
  transcript in \code{assembled}.
}
\description{
label assembled transcripts with gene names
}
\details{
chromosome labels in \code{gtf} and \code{assembled} should match.
  (i.e., you should provide the path to a gtf corrsponding to the same
  annotation you used when constructing \code{assembled})
}
\examples{
data(bg)
gtfPath = system.file('extdata', 'annot.gtf.gz', package='ballgown')
geneoverlaps = getGenes(gtfPath, structure(bg)$trans, UCSC=FALSE)
}
\author{
Alyssa Frazee, Andrew Jaffe
}
