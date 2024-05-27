<div align="center">
    <img src="http://146.56.237.198:3838/findGSEP/findGSEP_logo1205-2.png" width="50%" height=auto>
</div>

# findGSEP

<!-- badges: start -->
  
[![R-CMD-check](https://github.com/sperfu/findGSEP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sperfu/findGSEP/actions/workflows/R-CMD-check.yaml)

[![](https://www.r-pkg.org/badges/version/findGSEP?color=orange)](https://cran.r-project.org/package=findGSEP)

[![License: GPL2](https://img.shields.io/badge/license-GPL2-blue.svg)](https://cran.r-project.org/web/licenses/GPL2)
<!-- badges: end -->


Accurate estimating genome size is a crucial task in sequencing projects. Current methods often struggle with polyploidy or become inefficient when dealing with species that exceed a ploidy level of six. To address these challenges, we introduce findGSEP, an enhanced version of findGSE. findGSEP utilizes a segmented fitting approach to fit a normal distribution to polyploid species within a segmented framework. This ap-proach simplifies the process of single fitting while significantly expanding the range of ploidy levels it can handle. Moreover, findGSEP offers users interactive tools through both an open-source R application and a web application, facilitating reliable and precise estimation of genome size.  

# News 🌟

We have released our backend-server findGSEP and provide a CPU-based version of [findGSEP](http://146.56.237.198:3838/findGSEP/) online platform. Please check it out!!!


## Installation & Usage

### Instructions for running Jellyfish:

1. **Download and install jellyfish** from: [Jellyfish Release](http://www.genome.umd.edu/jellyfish.html#Release)

2. **Count kmers using jellyfish**:

    ```sh
    jellyfish count -C -m 21 -t 1 -s 5G *.fastq -o reads.mer
    ```
    
    *Note: Adjust the memory (-s) and threads (-t) parameters according to your server. This example uses 1 thread and 5GB of RAM. The kmer length (-m) may need to be scaled if you have low coverage or a high error rate. Always use 'canonical kmers' (-C).*

3. **Export the kmer count histogram**:

    ```sh
    jellyfish histo -h 3000000 -t 10 -o reads.histo reads.mer
    ```
    
    *Note: The thread count (-t) should be scaled according to your server.*

4. **Upload `reads.histo` to findGSEP**.


### Using KMC:

1. **Download and install KMC** from: [KMC GitHub](https://github.com/refresh-bio/KMC)

2. **Count kmers using KMC**:

    ```sh
    kmc -k21 -t1 -m5 -ci1 *.fastq reads_kmc tmp
    ```
    
    *Note: Adjust the memory (-m) and threads (-t) parameters according to your server. This example uses 1 thread and 5GB of RAM. The kmer length (-k) may need to be scaled if you have low coverage or a high error rate. The `-ci1` option ensures that kmers with a count of at least 1 are included.*

3. **Export the kmer count histogram**:

    ```sh
    kmc_tools transform reads_kmc histogram reads_kmc.histo
    ```
    
    *Note: This will create the histogram file `reads_kmc.histo`.*

4. **Upload `reads_kmc.histo` to findGSEP**.


## Instructions for installing findGSEP package

### Get the released version from CRAN:

``` r
install.packages("findGSEP")
```

### Or the development version from github:


1. Install `devtools`:

```bash
install.packages("devtools")
```

2. Install directly from GitHub:


```bash
devtools::install_github("sperfu/findGSEP")
```

#### Note: This package was developed using R version 4.2.0. To ensure the stability of the package, it is highly recommended that users install R version 4.2.0.


### Data

You can check our demo dataset at our [webserver](http://146.56.237.198:3838/findGSEP/) or [drive](https://drive.google.com/drive/folders/1fBuPRxi_J-oMpj6G2KokEUeB_8S8ahH6?usp=share_link) for complete data. We have provide precalculated histo file whose ploidy number ranging from tetraploid to octoploid.

### Usage:
```R
# Set options (optional):

options(warn = -1)

# Define input parameters:

path <- "histo_files"
samples <- "your_file.histo"
sizek <- 21
exp_hom <- 200
ploidy <- 4
output_dir <- "outfiles"
xlimit <- -1
ylimit <- -1
range_left <- exp_hom * 0.2
range_right <- exp_hom * 0.2

#Call the findGSEP function with specified parameters:

findGSEP(path, samples, sizek, exp_hom, ploidy, range_left, range_right, xlimit, ylimit, output_dir)

# For any questions, usage inquiries, or reporting potential bugs, please contact the author.
```

After running, You will find 'your_file.histo_hap_genome_size_est.pdf' in your output_dir folder, please give it a try!!!

## Parameter settings

You can reference to our paramenter setting for those species we used in our [webserver](http://146.56.237.198:3838/findGSEP/) or [demo dataset](https://drive.google.com/drive/folders/1fBuPRxi_J-oMpj6G2KokEUeB_8S8ahH6?usp=share_link).

| Species           | Expected Hom(Mb) | Ploidy number | Size k |
|-------------------|------------------|---------------|--------|
| Chinese sturgeon  | 100              | 8             | 21     |
| Strawberry        | 100              | 8             | 21     |
| Wheat             | 150              | 6             | 21     |
| Redwood           | 80               | 6             | 21     |
| Cotton            | 150              | 4             | 21     |
| Javanica          | 200              | 4             | 21     |
| Potato            | 180              | 4             | 21     |
| Floridensis       | 220              | 4             | 21     |
| Crayfish          | 35               | 3             | 21     |
| Enterolobii       | 130              | 3             | 21     |
| Incognita         | 200              | 3             | 21     |
| Seabass           | 80               | 2             | 21     |
| Bird              | 40               | 2             | 21     |
| Drosophila        | 50               | 2             | 21     |
| Pear              | 100              | 2             | 21     |
| Oyster            | 50               | 2             | 21     |


## Note:

If you enconter problem when installing devtools, especially for those packages below, please consider install them using conda install command:

```bash
conda install -c conda-forge r-gert
conda install -c conda-forge r-textshaping
conda install -c conda-forge r-ragg
conda install -c conda-forge r-pkgdown
```

If you enconter issues like:

1. could not find function "brewer.pal"
   
2.  could not find function "alpha"

Solutions:

```R
library(RColorBrewer)
library(ggplot2)
```



