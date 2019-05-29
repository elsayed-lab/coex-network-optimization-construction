# Installation

1. Create a conda environment:

```sh
conda update conda
conda create -n coex-nets -c bioconda -c defaults -c conda-forge r-base=3.5.1 r-devtools \
     r-rcppparallel r-rcpparmadillo r-rmarkdown gcc_linux-64 gfortran_linux-64
```

2. Initialize the renv environment

Launch R and run:

```
renv::init()
```

Re-initialize the environment (option 2).

Next, install Bioconductor dependencies using:

```r
BiocManager::install(c('Biobase', 'biomaRt', 'goseq', 'GO.db',  'KEGGREST', 'limma',
                       'impute', 'rtracklayer', 'preprocessCore', 'sva', 'WGCNA'))
```
