<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Simulating RNA-seq reads with Polyester}
-->

# The Polyester package for simulating RNA-seq reads

## Why use Polyester?

Polyester is an R package designed to simulate an RNA sequencing experiment. Given a set of annotated transcripts, polyester will simulate the steps of an RNA-seq experiment (fragmentation, reverse-complementing, and sequencing) and produce files containing simulated RNA-seq reads. Simulated reads can be analyzed using any of several downstream analysis tools. 

In particular, Polyester was designed to simulate a case/control experiment with biological replicates. Users are able to set differential transcript expression between cases and controls. This allows users to create datasets with known differential expression, which means they can the accuracy of statistical methods for differential expression detection.

Polyester was developed with several specific features in mind:  
* Simulation of differential expression at the transcript level
* Ability to set differential expression signal strength
* Simulation of small datasets, since large RNA-seq datasets can require lots of time and computing resources to analyze
* Generation of raw RNA-seq reads (as opposed to read alignments or transcript-level abundance estimates)
* Transparency: open-source code, 

## Prerequisites

Polyester depends on the `Biostrings` library from Bioconductor. You can install Biostrings by opening R and running:
```S
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```

We also recommend using R >= 3.0.0: because this vignette was written with [knitr](http://yihui.name/knitr/), it won't be compiled upon package installation with R versions < 3.0.0. (Support for non-Sweave vignettes was introduced in R 3.0.0). A vignette-less Polyester will likely work with older versions of R, but will not be officially supported.

## Simulating reads


## Future features
The following features will be implemented in a future release: 
* ability to set a seed, so running the simulation function twice will produce the same set of simulated reads (**near future**)
* option to simulate from GTF file + DNA sequence, instead of FASTA file of transcripts





We can use the JavaScript library [`DataTables`](http://www.datatables.net) to generate enhanced tables in HTML. In the example below, we create a table for the `mtcars` data:


```r
library(knitr)
kable(mtcars, "html", table.attr = "id=\"mtcars_table\"")
```

<table id="mtcars_table">
 <thead>
  <tr>
   <th>   </th>
   <th> mpg </th>
   <th> cyl </th>
   <th> disp </th>
   <th> hp </th>
   <th> drat </th>
   <th> wt </th>
   <th> qsec </th>
   <th> vs </th>
   <th> am </th>
   <th> gear </th>
   <th> carb </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> Mazda RX4 </td>
   <td> 21.0 </td>
   <td> 6 </td>
   <td> 160.0 </td>
   <td> 110 </td>
   <td> 3.90 </td>
   <td> 2.620 </td>
   <td> 16.46 </td>
   <td> 0 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Mazda RX4 Wag </td>
   <td> 21.0 </td>
   <td> 6 </td>
   <td> 160.0 </td>
   <td> 110 </td>
   <td> 3.90 </td>
   <td> 2.875 </td>
   <td> 17.02 </td>
   <td> 0 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Datsun 710 </td>
   <td> 22.8 </td>
   <td> 4 </td>
   <td> 108.0 </td>
   <td>  93 </td>
   <td> 3.85 </td>
   <td> 2.320 </td>
   <td> 18.61 </td>
   <td> 1 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Hornet 4 Drive </td>
   <td> 21.4 </td>
   <td> 6 </td>
   <td> 258.0 </td>
   <td> 110 </td>
   <td> 3.08 </td>
   <td> 3.215 </td>
   <td> 19.44 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Hornet Sportabout </td>
   <td> 18.7 </td>
   <td> 8 </td>
   <td> 360.0 </td>
   <td> 175 </td>
   <td> 3.15 </td>
   <td> 3.440 </td>
   <td> 17.02 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Valiant </td>
   <td> 18.1 </td>
   <td> 6 </td>
   <td> 225.0 </td>
   <td> 105 </td>
   <td> 2.76 </td>
   <td> 3.460 </td>
   <td> 20.22 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Duster 360 </td>
   <td> 14.3 </td>
   <td> 8 </td>
   <td> 360.0 </td>
   <td> 245 </td>
   <td> 3.21 </td>
   <td> 3.570 </td>
   <td> 15.84 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Merc 240D </td>
   <td> 24.4 </td>
   <td> 4 </td>
   <td> 146.7 </td>
   <td>  62 </td>
   <td> 3.69 </td>
   <td> 3.190 </td>
   <td> 20.00 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 4 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Merc 230 </td>
   <td> 22.8 </td>
   <td> 4 </td>
   <td> 140.8 </td>
   <td>  95 </td>
   <td> 3.92 </td>
   <td> 3.150 </td>
   <td> 22.90 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 4 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Merc 280 </td>
   <td> 19.2 </td>
   <td> 6 </td>
   <td> 167.6 </td>
   <td> 123 </td>
   <td> 3.92 </td>
   <td> 3.440 </td>
   <td> 18.30 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 4 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Merc 280C </td>
   <td> 17.8 </td>
   <td> 6 </td>
   <td> 167.6 </td>
   <td> 123 </td>
   <td> 3.92 </td>
   <td> 3.440 </td>
   <td> 18.90 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 4 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Merc 450SE </td>
   <td> 16.4 </td>
   <td> 8 </td>
   <td> 275.8 </td>
   <td> 180 </td>
   <td> 3.07 </td>
   <td> 4.070 </td>
   <td> 17.40 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 3 </td>
  </tr>
  <tr>
   <td> Merc 450SL </td>
   <td> 17.3 </td>
   <td> 8 </td>
   <td> 275.8 </td>
   <td> 180 </td>
   <td> 3.07 </td>
   <td> 3.730 </td>
   <td> 17.60 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 3 </td>
  </tr>
  <tr>
   <td> Merc 450SLC </td>
   <td> 15.2 </td>
   <td> 8 </td>
   <td> 275.8 </td>
   <td> 180 </td>
   <td> 3.07 </td>
   <td> 3.780 </td>
   <td> 18.00 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 3 </td>
  </tr>
  <tr>
   <td> Cadillac Fleetwood </td>
   <td> 10.4 </td>
   <td> 8 </td>
   <td> 472.0 </td>
   <td> 205 </td>
   <td> 2.93 </td>
   <td> 5.250 </td>
   <td> 17.98 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Lincoln Continental </td>
   <td> 10.4 </td>
   <td> 8 </td>
   <td> 460.0 </td>
   <td> 215 </td>
   <td> 3.00 </td>
   <td> 5.424 </td>
   <td> 17.82 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Chrysler Imperial </td>
   <td> 14.7 </td>
   <td> 8 </td>
   <td> 440.0 </td>
   <td> 230 </td>
   <td> 3.23 </td>
   <td> 5.345 </td>
   <td> 17.42 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Fiat 128 </td>
   <td> 32.4 </td>
   <td> 4 </td>
   <td>  78.7 </td>
   <td>  66 </td>
   <td> 4.08 </td>
   <td> 2.200 </td>
   <td> 19.47 </td>
   <td> 1 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Honda Civic </td>
   <td> 30.4 </td>
   <td> 4 </td>
   <td>  75.7 </td>
   <td>  52 </td>
   <td> 4.93 </td>
   <td> 1.615 </td>
   <td> 18.52 </td>
   <td> 1 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Toyota Corolla </td>
   <td> 33.9 </td>
   <td> 4 </td>
   <td>  71.1 </td>
   <td>  65 </td>
   <td> 4.22 </td>
   <td> 1.835 </td>
   <td> 19.90 </td>
   <td> 1 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Toyota Corona </td>
   <td> 21.5 </td>
   <td> 4 </td>
   <td> 120.1 </td>
   <td>  97 </td>
   <td> 3.70 </td>
   <td> 2.465 </td>
   <td> 20.01 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Dodge Challenger </td>
   <td> 15.5 </td>
   <td> 8 </td>
   <td> 318.0 </td>
   <td> 150 </td>
   <td> 2.76 </td>
   <td> 3.520 </td>
   <td> 16.87 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> AMC Javelin </td>
   <td> 15.2 </td>
   <td> 8 </td>
   <td> 304.0 </td>
   <td> 150 </td>
   <td> 3.15 </td>
   <td> 3.435 </td>
   <td> 17.30 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Camaro Z28 </td>
   <td> 13.3 </td>
   <td> 8 </td>
   <td> 350.0 </td>
   <td> 245 </td>
   <td> 3.73 </td>
   <td> 3.840 </td>
   <td> 15.41 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Pontiac Firebird </td>
   <td> 19.2 </td>
   <td> 8 </td>
   <td> 400.0 </td>
   <td> 175 </td>
   <td> 3.08 </td>
   <td> 3.845 </td>
   <td> 17.05 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Fiat X1-9 </td>
   <td> 27.3 </td>
   <td> 4 </td>
   <td>  79.0 </td>
   <td>  66 </td>
   <td> 4.08 </td>
   <td> 1.935 </td>
   <td> 18.90 </td>
   <td> 1 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Porsche 914-2 </td>
   <td> 26.0 </td>
   <td> 4 </td>
   <td> 120.3 </td>
   <td>  91 </td>
   <td> 4.43 </td>
   <td> 2.140 </td>
   <td> 16.70 </td>
   <td> 0 </td>
   <td> 1 </td>
   <td> 5 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Lotus Europa </td>
   <td> 30.4 </td>
   <td> 4 </td>
   <td>  95.1 </td>
   <td> 113 </td>
   <td> 3.77 </td>
   <td> 1.513 </td>
   <td> 16.90 </td>
   <td> 1 </td>
   <td> 1 </td>
   <td> 5 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Ford Pantera L </td>
   <td> 15.8 </td>
   <td> 8 </td>
   <td> 351.0 </td>
   <td> 264 </td>
   <td> 4.22 </td>
   <td> 3.170 </td>
   <td> 14.50 </td>
   <td> 0 </td>
   <td> 1 </td>
   <td> 5 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Ferrari Dino </td>
   <td> 19.7 </td>
   <td> 6 </td>
   <td> 145.0 </td>
   <td> 175 </td>
   <td> 3.62 </td>
   <td> 2.770 </td>
   <td> 15.50 </td>
   <td> 0 </td>
   <td> 1 </td>
   <td> 5 </td>
   <td> 6 </td>
  </tr>
  <tr>
   <td> Maserati Bora </td>
   <td> 15.0 </td>
   <td> 8 </td>
   <td> 301.0 </td>
   <td> 335 </td>
   <td> 3.54 </td>
   <td> 3.570 </td>
   <td> 14.60 </td>
   <td> 0 </td>
   <td> 1 </td>
   <td> 5 </td>
   <td> 8 </td>
  </tr>
  <tr>
   <td> Volvo 142E </td>
   <td> 21.4 </td>
   <td> 4 </td>
   <td> 121.0 </td>
   <td> 109 </td>
   <td> 4.11 </td>
   <td> 2.780 </td>
   <td> 18.60 </td>
   <td> 1 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 2 </td>
  </tr>
</tbody>
</table>


Note we assigned an `id` to the table, and next we use the `DataTables` library to tweak the table.

```js
<script type="text/javascript" charset="utf-8">
    $(document).ready(function() {
        $('#mtcars_table').dataTable();
    } );
</script>
```

<script type="text/javascript" charset="utf-8">
  $(document).ready(function() {
        $('#mtcars_table').dataTable();
    } );
</script>

Since this is a Markdown vignette, we need to add the JavaScript libraries as well as some additional CSS files to the HTML header, and this can be done via:


```r
options(markdown.HTML.header = system.file("misc", "datatables.html", package = "knitr"))
```


Two files `vignette.css` and `datatables.html` under the `misc` directory of **knitr** were added to the HTML output.

By comparison, below is an ordinary table:


```r
kable(head(mtcars), "html")
```

<table>
 <thead>
  <tr>
   <th>   </th>
   <th> mpg </th>
   <th> cyl </th>
   <th> disp </th>
   <th> hp </th>
   <th> drat </th>
   <th> wt </th>
   <th> qsec </th>
   <th> vs </th>
   <th> am </th>
   <th> gear </th>
   <th> carb </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> Mazda RX4 </td>
   <td> 21.0 </td>
   <td> 6 </td>
   <td> 160 </td>
   <td> 110 </td>
   <td> 3.90 </td>
   <td> 2.620 </td>
   <td> 16.46 </td>
   <td> 0 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Mazda RX4 Wag </td>
   <td> 21.0 </td>
   <td> 6 </td>
   <td> 160 </td>
   <td> 110 </td>
   <td> 3.90 </td>
   <td> 2.875 </td>
   <td> 17.02 </td>
   <td> 0 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 4 </td>
  </tr>
  <tr>
   <td> Datsun 710 </td>
   <td> 22.8 </td>
   <td> 4 </td>
   <td> 108 </td>
   <td>  93 </td>
   <td> 3.85 </td>
   <td> 2.320 </td>
   <td> 18.61 </td>
   <td> 1 </td>
   <td> 1 </td>
   <td> 4 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Hornet 4 Drive </td>
   <td> 21.4 </td>
   <td> 6 </td>
   <td> 258 </td>
   <td> 110 </td>
   <td> 3.08 </td>
   <td> 3.215 </td>
   <td> 19.44 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 1 </td>
  </tr>
  <tr>
   <td> Hornet Sportabout </td>
   <td> 18.7 </td>
   <td> 8 </td>
   <td> 360 </td>
   <td> 175 </td>
   <td> 3.15 </td>
   <td> 3.440 </td>
   <td> 17.02 </td>
   <td> 0 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 2 </td>
  </tr>
  <tr>
   <td> Valiant </td>
   <td> 18.1 </td>
   <td> 6 </td>
   <td> 225 </td>
   <td> 105 </td>
   <td> 2.76 </td>
   <td> 3.460 </td>
   <td> 20.22 </td>
   <td> 1 </td>
   <td> 0 </td>
   <td> 3 </td>
   <td> 1 </td>
  </tr>
</tbody>
</table>

