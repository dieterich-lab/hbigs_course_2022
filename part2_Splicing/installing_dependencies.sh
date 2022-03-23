#!/usr/bin/env bash


# python libs
pip install snakemake --user
pip install jupyter --user

# R packages
Rscript -e 'install.packages(c("BiocManager", "knitr", "devtools",  "R.utils", "Hmisc", "dplyr", "stringr", "IRkernel"), lib="~/R/x86_64-pc-linux-gnu-library/3.6", repos="https://cran.uni-muenster.de", INSTALL_opts = c("--no-lock"))'

# IRkernel
Rscript -e 'IRkernel::installspec()'

# DRIMSeq
Rscript -e 'BiocManager::install(c("DRIMSeq", "rtracklayer", "tximport"), INSTALL_opts = c("--no-lock"))'