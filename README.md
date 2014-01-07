# Ballgown
R package for downstream analysis of transcript assemblies. This is the bleeding-edge, development version.  Please do submit issues, pull requests, bug fixes, usability improvements, etc.

## (0) installation
After installing the `devtools` package from CRAN and the `GenomicRanges` package from Bioconductor, run:

```S
library(devtools)
install_github('ballgown', 'alyssafrazee')
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

```S
library(ballgown)
dirs <- c('sample01_output', 'sample02_output')
awesome_bg <- ballgown(dirs = dirs)
save(awesome_bg, file = 'awesome_bg.rda')
```

If you have a big experiment, loading the data might require a lot of memory and a lot of time, so it might be best to do this as a non-interactive batch job.  You only have to do this once, though.  The resulting `rda` file is usually only a few Gb on disk, even for large experiments, and usually only takes a reasonable amount of memory to work with.  (It's the creation of the object that's the memory-hog).  After initial creation of the object, whenever you want to work with this assembly, you can just fire up R and load it in:

```S
library(ballgown)
load('awesome_bg.rda')
awesome_bg
# ballgown instance with 109283 assembled transcripts
```
## (3) accessing assembly data
A `ballgown` object has three main components: `structure`, `data`, and `indexes`.  The `structure` component,  depending heavily on the `GenomicRanges` Bioconductor package, specifies the strucutre (i.e., genomic locations and relationships between exons, introns, and transcripts) of your assembly.  It's convenient to represent exons and introns as intervals and to represent transcripts as a set of intervals (exons), so the assembled exons and introns are available as `GRanges` objects, and the assembled transcripts are available as a `GRangesList` object.  This means that useful range operations, such as `findOverlaps` and `reduce`, are readily available for your assembled features. You can easily extract these objects from the main `ballgown` object or examine them directly:

```S
structure(awesome_bg)$exon
structure(awesome_bg)$intron
transcript_struct <- structure(awesome_bg)$trans
```
The second main component of a `ballgown` object is `data`, i.e., tables containing expression data for the genomic features.  These tables are very similar to the `*_data.ctab` tables described in section (1).  In general, you can extract expression data using the syntax `*expr(ballgown_object_name, <EXPRESSION_MEASUREMENT>)`, where `*` is either e for exon, i for intron, t for transcript, or g for gene, and <EXPRESSION MEASUREMENT> is an expression-measurement column name from the appropriate `.ctab` file.  Gene-level measurements are calculated by appropriately aggregating the transcript-level measurements for that gene.  All of the following are valid ways to extract expression data from the `awesome_bg` ballgown object:

```S
transcript_fpkm <- texpr(awesome_bg, 'FPKM')
transcript_cov <- texpr(awesome_bg, 'cov')
whole_tx_table <- texpr(awesome_bg)
exon_mcov <- eexpr(awesome_bg, 'mcov')
junction_cov <- iexpr(awesome_bg, 'cov')
whole_intron_table <- iexpr(awesome_bg)
gene_expression <- gexpr(awesome_bg)
```
Calculating the gene-level expression measurements can be slow for large experiments, so you may want to run the `gexpr` call as a batch job and save the result as an rda file for later use. 

Finally, the `indexes` component of the ballgown object connects the various pieces of the assembly and provides other information about your experiment.  Importantly, there is a slot called `pData` that holds a data frame of phenotype information for the samples.  Usually you have to create this manually.  **Make sure of two things: (a) the column of this data frame that identifies the samples is called `dirname` and (b) that column is ordered the same way as the tables in the `data` component.**  You can check the order by running something like `names(texpr(awesome_bg))`, or you can ensure the sample directories are in alphabetical order (i.e., in the order they'd appear if you ran an `ls` on the root directory). 

You can assign `pData` in several different ways, two of which are shown below.  It's probably easiest to assign in as you load the data, so you don't have to re-save the object, but if you forget, no biggie.

```S
## creating pData as you create the object:
dirs <- c('sample01_output', 'sample02_output')
pData <- data.frame(dirname = dirs, population = c('normal', 'cancer'))
awesome_bg <- ballgown(dirs = dirs, pData = pData)

## creating pData after the object has been created:
load('awesome_bg.rda')
dirs <- c('sample01_output', 'sample02_output')
pData(awesome_bg) <- data.frame(dirname = dirs, population = c('normal', 'cancer'))
save(awesome_bg, file = 'awesome_bg.rda') #re-save the object with associated pData
```
The other components of `indexes` are the `e2t` and `i2t` tables described in section (1), as well as a `t2g` table denoting which transcripts belong to which genes.  There is also a `bamfiles` component, designed to hold the file paths to the read alignment files for each sample.  Again, make sure this is in the same order as the `dirname` column of `pData`.  The `bamfiles` component isn't currently used by any ballgown functions, but we imagine it could be useful for fans of `RSamtools` or similar packages.  Here are some examples of how to extract `indexes` components from ballgown objects:

```S
exon_transcript_table <- indexes(awesome_bg)$e2t
transcript_gene_table <- indexes(awesome_bg)$t2g
alignment_files <- indexes(awesome_bg)$bamfiles
phenotype_table <- pData(awesome_bg)
```



## (4) visualizing transcript structure

You can see what the assembled transcripts look like for each gene using the `plotTranscripts` function.  Transcripts or exons can be colored by expression level.  Features can also be colored by mean expression level for a population or for the entire experiment.  The possibilities are endless!  

```S
plotTranscripts(gene = 'XLOC_000043', gown = awesome_bg, samp = 'FPKM.sample01_output', colorby = 'transcript', main = 'transcripts from gene XLOC_000043: sample 1, FPKM')
```

![example plot](https://raw.github.com/alyssafrazee/ballgown/master/explot.png)


## (5) differential expression analysis

You can choose from a wide selection of simple, fast statistical methods for testing whether transcripts are differentially expressed between experimental conditions or across time.  The main way to do this in the ballgown package is with the flexible `stattest` function.   This function currently supports two-group (e.g., case/control) comparisons, multi-group comparisons, and timecourse differential expression testing.  For multi-group comparisons, a significant result would indicate that the feature is differentially expressed in at least one of the groups, while for timecourse comparisons, a significant result would indicate that the feature has an expression profile that varies over time (versus an expression profile that is flat, or contstant, over time).

The default statistical test in ballgown is a parametric F-test comparing nested linear models.  This is described in detail in our paper (coming soon!), but briefly, two models are fit for each transcript: one including the covariate of interest (e.g., case/control status or time) and one not including that covariate.  An [F statistic](http://en.wikipedia.org/wiki/F-test#Regression_problems) is calculated using the fits of the two models, and a corresponding p-value is obtained.  A significant p-value means the model including the covariate of interest fits significantly better than the model without that covariate, which indicates differential expression.  We adjust for multiple testing by reporting q-values for each transcript in addition to p-values: reporting features with, say, q < 0.05 means the false discovery rate should be controlled at about 5%. 

Here is an example of how to use `stattest`:
```S
stat_results <- stattest(awesome_bg, feature = 'transcript', meas = 'FPKM', covariate = 'population')
```
The resulting object, `stat_results`, is a data frame containing the feature tested, feature ids, and corresponding p- and q-values.

At minimum, you need to provide the name of the ballgown object, which type of feature you want to test ( gene, transcript, exon, or intron), the expression measurement you want to use (FPKM, cov, rcount, etc.), and the covariate of interest, which must be the name of one of the columns of the `pData` component of your ballgown object.  This variable is automatically converted to a factor during model fitting.

### timecourse experiments
For timecourse experiments, we fit a smooth curve to time using [natural splines](http://en.wikipedia.org/wiki/Spline_interpolation).  The model including these spline terms is compared to a model without any spline terms for the F-test.  A typical call to `stattest` for a timecourse experiment would look like:
```S
timecourse_results <- stattest(awesome_bg, feature = 'transcript', meas = 'FPKM', covariate = 'time', timecourse = TRUE)
```
In the timecourse scenario, `covariate` refers to a numeric time variable in `pData`, and you must specify `timecourse = TRUE`.  

### adjusting for confounders
You can adjust for any or all variables in `pData` when testing for differential expression.  We automatically adjust for library size using the sum of all expression measurements below the 75th percentile for that sample, but we hope to allow users to choose their own library-size adjustment in the near future.  If you would like to adjust for other confounders, say sex and age, just provide those confounders as the `adjustvars` argument to `stattest`, e.g., `adjustvars = c('age', 'sex')`.

### using alternative statistical methods
Our statistical methods for differential expression testing are straightforward and accurate, but if you would rather use another tool, that's cool too!  Ballgown's data structures make it easy to use table-based packages like [limma](http://www.bioconductor.org/packages/2.13/bioc/html/limma.html), [limma Voom](http://www.statsci.org/smyth/pubs/VoomPreprint.pdf), [DESeq](http://www.bioconductor.org/packages/release/bioc/html/DESeq.html), [DEXSeq](http://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html), or [EdgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) for differential expression analysis.  A transcript-by-sample expression table can be easily created with, say, `my_table <- texpr(awesome_bg, 'cov')`, and `my_table` can be used as the input for these or other packages.

## that's it!
Functions are completely documented, so help is available with `help(function)` or `?function`.  And if you find a bug or other problem, we definitely want to know about it!

## footnotes
1.  email me (acfrazee@gmail.com) for a copy of the `tablemaker` binary.  It should be public shortly, but for now, it's an "upon request" situation.