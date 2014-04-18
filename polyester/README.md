# The Polyester package for simulating RNA-seq reads

## Why use Polyester?

Polyester is an R package designed to simulate an RNA sequencing experiment. Given a set of annotated transcripts, polyester will simulate the steps of an RNA-seq experiment (fragmentation, reverse-complementing, and sequencing) and produce files containing simulated RNA-seq reads. Simulated reads can be analyzed using any of several downstream analysis tools. 

In particular, Polyester was designed to simulate a case/control experiment with biological replicates. Users are able to set differential transcript expression between cases and controls. This allows users to create datasets with known differential expression, which means they can the accuracy of statistical methods for differential expression detection.

Polyester was developed with several specific features in mind:  
* Simulation of differential expression at the transcript level
* Ability to set differential expression signal strength
* Simulation of small datasets, since large RNA-seq datasets can require lots of time and computing resources to analyze
* Generation of raw RNA-seq reads (as opposed to read alignments or transcript-level abundance estimates)
* Transparency/open-source code

## Prerequisites

Polyester depends on the `Biostrings` and `IRanges` libraries from Bioconductor. You can install these packagess by opening R and running:
```S
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("IRanges")
```

We also recommend using R >= 3.0.0: because this vignette was written with [knitr](http://yihui.name/knitr/), it won't be compiled upon package installation with R versions < 3.0.0. (Support for non-Sweave vignettes was introduced in R 3.0.0). A vignette-less Polyester will likely work with older versions of R, but will not be officially supported.

Finally, you will need a reference FASTA file containing names and sequences of transcripts from which reads should be simulated. Known transcripts from human chromosome 22 (hg19 build) are available in the `data` subdirectory of this package. 

## Simulating reads

Simulating an RNA-seq experiment with Polyester requires just one function call. You can choose either `simulate_experiment()` or `simulate_experiment_countmat()`. See the function-specific documentation for examples on using each one. The ideas behind the two approaches are:

### approach 1: built-in negative binomial model (two-group experiment)
The `simulate_experiment` function draws the number of reads to simulate from each transcript from a negative binomial distribution. For this function, you need to specify:
* `num_reps`: Number of biological replicates per experimental group (default: 10; can specify different numbers of replicates in the groups)
* `fold_changes`: A fold change for each transcript. This fold change represents the multiplicative change in the mean number of reads generated from each transcript, between the two experimental groups.
* `reads_per_transcript`: The baseline mean number of reads for each transcript. 
    - Fold changes compare the mean number of reads in group 1 to group 2. So a fold change of 0.5 means group 2's baseline mean number of reads for this transcript is twice that of group 1.
    - Long transcripts usually produce more reads in RNA-seq experiments than short ones, so you may want to specify `reads_per_transcript` as a function of transcript length
    - Default is 300 (regardless of transcript length).
* `dispersion_param`: controls the per-transcript mean/variance relationship. In the negative binomial distribution, the mean/variance relationship is: ```mean = mean + (mean^2) / size```, where "size" is the dispersion parameter. You can specify the dispersion parameter for each transcript. By default, size is defined as 1/3 of the transcript's mean.

### approach 2: build your own expression model
The `simulate_experiment_readmat` function takes a count matrix as an argunent. Each row of this matrix represents a transcript, and each column represents a sample in the experiment. Entry `i,j` of the matrix specifies how many reads should be sampled from transcript `i` for sample `j`, allowing you to precisely and flexibly define the (differential) transcript expression structure for the experiment.

### other simulation parameters that can be set:
For both `simulate_experiment` and `simulate_experiment_countmat`, you can change these parameters:
* `fraglen`: Mean fragment length (default 250)
* `fragsd`: Standard devation of fragment lengths (default 25)
* `readlen`: Read length (default 100)
* `error_rate`: Sequencing error rate: probability that the sequencer records the wrong nucleotide at any given base (default 0.005, uniform error model assumed)
* `paired`: Whether the reads should be paired-end (default TRUE)

[This review paper](http://genomebiology.com/2010/11/12/220) (Oshlack, Robinson, and Young, _Genome Biology_ 2010, open access) provides a good overview of the RNA sequencing process, and might be particularly useful for understanding where some of these simulation parameters come into play.

If you'd like to explore specific steps in the sequencing process (fragmentation, reverse-complementing, error-adding), the functions called within `simulate_experiment` are also available and individually documented in Polyester.

### output
A call to `simulate_experiment` or `simulate_experiment_countmat` will write FASTA files to the directory specified by the `outdir` argument. Reads in the FASTA file will be labeled with the transcript from which they were simulated.

If `paired` is true, you'll get two FASTA files per biological replicate (left mates are designated by the suffix `_1.fasta`; right mates by `_2.fasta`). If single-end reads are generated (`paired=FALSE`) you'll get one FASTA file per replicate. 

Files will be named `sample_01` through `sample_N` where `N` is the total number of replicates. The first `num_reps` (or `num_reps[1]`) samples belong to the same group in the two-group experiment scenario. 

In `simulate_experiment`, by default, a table called `sim_info.txt` is written to `outdir`, which will contain transcript IDs, fold changes, and whether or not that transcript was set to be differentially expressed. This file could be useful for downstream analysis. If the transcript names in the FASTA file cause problems down the line (e.g., a dangling single quote from a `5'-end` label), you can specify your own transcript names with the `transcriptid` argument. You will need to keep track of this information separately if you use `simulate_experiment_countmat.`

## Future features
The following features will be implemented in a future release: 
* ability to set a seed, so running the simulation function twice will produce the same set of simulated reads (**near future**)
* option to simulate from GTF file + DNA sequence, instead of FASTA file of transcripts

## Bug reports
Report bugs as issues on our [GitHub repository](https://github.com/alyssafrazee/ballgown/tree/master/polyester). 
