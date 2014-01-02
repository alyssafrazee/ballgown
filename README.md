# Ballgown
R package for downstream analysis of transcript assemblies. This is the bleeding-edge, development version.  Please do submit issues, pull requests, bug fixes, usability improvements, etc.

## (0) installation
After installing the `devtools` package from CRAN and the `GenomicRanges` package from Bioconductor, run:

```
library(devtools)
install_github('ballgown', 'alyssafrazee', subdir = 'ballgown')
```

## (1) preprocessing: preparing transcript assemblies for analysis
This package was created as a statistical backend for existing transcript assembly software.  In its current state, `ballgown` only provides support for assemblies created with [Cufflinks](http://cufflinks.cbcb.umd.edu/index.html), though we hope to be able to process output from other assemblers in the near future.  

Users need to run the `tablemaker` binary<sup>1</sup> to get their assembly output organized into a format that `ballgown` can load.  This program needs to be run on each RNA-seq sample in your experiment.  It requires your merged assembly, in `gtf` format, and the read alignments for each sample, in `bam` format.  You use `tablemaker` as follows:

`tablemaker -p 4 -q -W -G merged.gtf -o sample01_output read_alignments.bam`

* `-p` denotes how many threads to use (the program can take a few hours to run, but can be parallelized)
* The `-q` can be removed for more verbose output messages  
* `-W` and `-G merged.gtf` are required.  The `-W` is just telling the program to run in tablemaker mode, and the `-G` argument points to the merged assembly file, which gives the assembled transcripts' structures
* The argument to `-o` is the desired output directory for the sample (each sample should have its own output directory)
* The read alignment file is the last argument.  If reads were aligned with TopHat, this is usually some variant of `accepted_hits.bam`

The output is 5 files, written to the specified output directory:

* `e_data.ctab`: exon-level expression measurements.  One row per exon.  Columns are `e_id` (numeric exon id), `chr`, `strand`, `start`, `end` (genomic location of the exon), and the following expression measurements for each sample:
	* `rcount`:  reads overlapping the exon 
	* `ucount`: uniquely mapped reads overlapping the exon 
	* `mrcount`: multi-map-corrected number of reads overlapping the exon
	* `cov` average per-base read coverage 
	* `cov_sd`: standard deviation of per-base read coverage
	* `mcov`: multi-map-corrected average per-base read coverage
	* `mcov_sd`: standard deviation of multi-map-corrected per-base coverage
* `i_data.ctab`: intron- (i.e., junction-) level expression measurements.  One row per intron.  Columns are `i_id` (numeric intron id), `chr`, `strand`, `start`, `end` (genomic location of the intron), and the following expression measurements for each sample: 
	* `rcount`: number of reads supporting the intron
	* `ucount`: number of uniquely mapped reads supporting the intron
	* `mrcount`: multi-map-corrected number of reads supporting the intron
* `t_data.ctab`: transcript-level expression measurements.  One row per transcript.  Columns are:
	* `t_id`: numeric transcript id
	* `chr`, `strand`, `start`, `end`: genomic location of the transcript
	* `t_name`: Cufflinks-generated transcript id
	* `num_exons`: number of exons comprising the transcript
	* `length`: transcript length, including both exons and introns
	* `gene_id`: gene the transcript belongs to
	* `gene_name`: HUGO gene name for the transcript, if known
	* `cov`: per-base coverage for the transcript (available for each sample)
	* `FPKM`: Cufflinks-estimated FPKM for the transcript (available for each sample)
* `e2t.ctab`: table with two columns, `e_id` and `t_id`, denoting which exons belong to which transcripts.  These ids match the ids in the `e_data` and `t_data` tables.
* `i2t.ctab`: table with two columns, `i_id` and `t_id`, denoting which introns belong to which transcripts.  These ids match the ids in the `i_data` and `t_data` tables.


## (2) loading assemblies into R
Let's assume you've run `tablemaker` on all your samples and that each sample's ballgown output directory is a subfolder of the same root directory.  (This is pretty much required).  So we might have a folder structure like:

```
awesome_experiment/
├── sample01_output/
    ├── e2t.ctab
    ├── e_data.ctab
    ├── i2t.ctab
    ├── i_data.ctab
    └── t_data.ctab
├── sample02_output/
    ├── e2t.ctab
    ├── e_data.ctab
    ├── i2t.ctab
    ├── i_data.ctab
    └── t_data.ctab
```
Then, in your R session, assuming your working directory is `awesome_experiment`, you can run:

```
library(ballgown)
dirs <- c('sample01_output', 'sample02_output')
awesome_bg <- ballgown(dirs = dirs)
save(awesome_bg, file = 'awesome_bg.rda')
```

If you have a big experiment, loading the data might require a lot of memory and a lot of time, so it might be best to do this as a non-interactive batch job.  You only have to do this once, though.  The resulting `rda` file is usually only a few Gb on disk, even for large experiments, and usually only takes a reasonable amount of memory to work with.  (It's the creation of the object that's the memory-hog).  After initial creation of the object, whenever you want to work with this assembly, you can just fire up R and load it in:

```R
library(ballgown)
load('awesome_bg.rda')
awesome_bg
```

## visualizing transcript structure

There are some great plotting functions here!

## differential expression analysis

You can choose from a wide selection of simple, fast statistical methods for testing whether transcripts are differentially expressed between experimental conditions or across time.

## footnotes
1.  email me (acfrazee@gmail.com) for a copy of the `tablemaker` binary.  It should be public shortly, but for now, it's an "upon request" situation.