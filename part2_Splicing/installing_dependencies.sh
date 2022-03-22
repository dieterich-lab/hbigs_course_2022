#!/usr/bin/env bash


# python libs
pip install snakemake --user
pip install jupyter --user

# R packages
Rscript -e 'install.packages(c("BiocManager", "knitr", "devtools", "rstan",  "rstantools", "R.utils", "Hmisc", "foreach", "dplyr", "doMC", "optparse", "shinyjs", "intervals", "shinycssloaders", "IRkernel"), lib="~/R/x86_64-pc-linux-gnu-library/3.6", repos="https://cran.uni-muenster.de", INSTALL_opts = c("--no-lock"))'

# IRkernel
Rscript -e 'IRkernel::installspec()'

# DRIMSeq
Rscript -e 'BiocManager::install(c("DRIMSeq", "rtracklayer", "tximport"), INSTALL_opts = c("--no-lock"))'

# LeafCutter

# First, R stan
mkdir ~/.R/
echo 'CXX14 = g++ -std=c++1y
CXX14FLAGS = -O3 -Wno-unused-variable -Wno-unused-function -fPIC' >>  ~/.R/Makevars

Rscript -e 'install.packages("rstan", type = "source",  lib="~/R/x86_64-pc-linux-gnu-library/3.6", repos="https://cran.uni-muenster.de", INSTALL_opts = c("--no-lock"))'

git clone https://github.com/davidaknowles/leafcutter

cd leafcutter/leafcutter
R CMD INSTALL --build . --no-lock -l ~/R/x86_64-pc-linux-gnu-library/3.6/
echo 'export PATH="$HOME/leafcutter/scripts/:$PATH"' >> .bashrc 

cd part_2/

wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz 
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i salmon_index  -p 10 2> salmon_index.log


