\docType{class}
\name{ballgown-class}
\alias{Ballgown}
\alias{ballgown-class}
\title{Ballgown}
\description{
S4 class for storing and manipulating expression data from
assembled transcriptomes
}
\section{Slots}{

\describe{
\item{\code{data}}{tables containing expression data for genomic features (introns, exons, transcripts)}

\item{\code{structure}}{genomic locations of features and their relationships to one another}

\item{\code{indexes}}{tables connecting components of the assembly and providing other experimental information (e.g., phenotype data and locations of read alignment files)}

\item{\code{dirs}}{directories holding data created by \code{tablemaker}}

\item{\code{mergedDate}}{date the ballgown object was created}
}}
\author{
Alyssa Frazee, Leonardo Collado Torres, Jeff Leek
}

