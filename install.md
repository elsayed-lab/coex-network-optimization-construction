# Installation

1. Create a conda environment:

```sh
conda update conda
conda create -n coex-nets -c bioconda -c defaults -c conda-forge \
   r-base=3.6.0 
```

```sh
conda update conda
conda create -n coex-nets -c bioconda -c defaults -c conda-forge \
     gcc_linux-64 gfortran_linux-64 pandoc pandoc-citeproc \
     r-base=3.6.0 r-devtools r-rcppparallel r-rcpparmadillo r-rmarkdown \
     r-biocmanager r-wgcna r-foreach r-doparallel r-tidyverse r-reshape2 \
     r-rcurl r-gridextra r-annotables \
     bioconductor-biobase bioconductor-keggrest bioconductor-goseq \
     bioconductor-go.db bioconductor-preprocesscore bioconductor-limma \
     bioconductor-biomart bioconductor-impute bioconductor-sva \
     bioconductor-rtracklayer \
```

2. Initialize the renv environment

Launch R and re-initialize the renv environment using:

```
renv::init()
```

(choose option 2, "Re-initialize the project, discovering and installing R package
dependencies as required").

Next, install Bioconductor dependencies using:

```r
BiocManager::install(c('Biobase', 'biomaRt', 'goseq', 'GO.db',  'KEGGREST', 'limma',
                       'impute', 'rtracklayer', 'preprocessCore', 'sva', 'WGCNA'))
```

Finally, install the necessary annotation packages using:

```r
renv::install('elsayed-lab/org.LmjF.tritryp.db')
renv::install('elsayed-lab/org.TcCLB.esmer.tritryp.db')
renv::install('elsayed-lab/TxDb.TcruziCLBrenerEsmer.tritryp32.genes')
renv::install('elsayed-lab/TxDb.LmajorFriedlin.tritryp32.genes') 
renv::install('elsayed-lab/Leishmania.major.Friedlin')
renv::install('elsayed-lab/Trypanosoma.cruzi.CLBrener.Esmeraldo') 
```
