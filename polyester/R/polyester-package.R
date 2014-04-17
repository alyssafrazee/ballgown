#' Polyester: simulating RNA-seq reads including differential expression
#'
#' Polyester is an R package designed to simulate an RNA sequencing experiment. Given a set of 
#' annotated transcripts, polyester will simulate the steps of an RNA-seq experiment 
#' (fragmentation, reverse-complementing, and sequencing) and produce files containing 
#' simulated RNA-seq reads. Simulated reads can be analyzed using any of several downstream 
#' analysis tools.
#'
#' A single function call produces RNA-seq reads in FASTA format from a case/control experiment
#' including biological replicates. Differential expression between cases and controls can be set
#' by the user, facilitating comparisons of statistical differential expression methods for RNA-seq
#' data. See detailed documentation for \code{\link{simulate_experiment}} and 
#' \code{\link{simulate_experiment_countmat}}. 
#' 
#' A vignette is available by typing \code{browseVignettes("polyester")} in the R prompt.
#' 
#' @references Alyssa C Frazee, Geo Pertea, Andrew E Jaffe, Ben Langmead, Steven L Salzberg, 
#' Jeffrey T Leek (2014). Flexible isoform-level differential expression analysis with Ballgown.
#' BioRxiv preprint: \url{http://biorxiv.org/content/early/2014/03/30/003665}. 
#' @import Biostrings
#' @docType package
#' @name polyester
NULL



