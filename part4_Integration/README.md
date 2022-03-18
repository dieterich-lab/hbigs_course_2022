# HBIGS course 2021/22: Introduction to computational RNA biology
## Part 4 - Proteomics and Integrative Analysis

## Description

This directory contains material for the processing and analysis of the *Proteomics* data (**1_Proteomics_Analysis**) as well as *Integrative Multi-Omics* (**2_Integrative_Analysis**).

The analysis steps have been applied over the _EGF-driven protein synthesis_ case-study data from [D.A. Rothenberg et al. A Proteomics Approach to Profiling the Temporal Translational Response to Stress and Growth. iScience. 2018; 9:367-381](https://www.sciencedirect.com/science/article/pii/S2589004218301949?via%3Dihub).

Analysis should be performed in the following order:

### 1.  Proteomics
Here are provided scripts for running a standard pipeline for Processing and Analysis of Proteomics Data. The pipeline shows how we can transform/normalize the raw intensities, filter for reads quantified in a low number of experiments, imputing missing values, identify and remove batch effects, interpret the data, identify significantly abundant reads, etc..

### 2.  Integration
Here are provided the scripts for Integrated Multi-Omics analysis of the Proteomics and Transcriptomics data. The analyses provided include correlation between the two data types, clustering analysis, pathway enrichment, and network analyses.

## Dependencies

Before running the scripts, please make sure to have installed the following R-Package dependencies:

#### CRAN
**1.**  [readr](https://cran.r-project.org/web/packages/readr/index.html)
**2.**  [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
**3.**  [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)
**4.**  [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
**5.**  [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html#installation)
**6.**  [knitr](https://www.r-project.org/nosvn/pandoc/knitr.html)
**7.**  [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
**8.**  [igraph](https://cran.r-project.org/web/packages/igraph/index.html)
**9.**  [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
**10.** [M2SMF](https://cran.r-project.org/web/packages/M2SMF/index.html)
**11.** [SNFtool](https://cran.r-project.org/web/packages/SNFtool/index.html)
**12.** [GSA](https://cran.r-project.org/web/packages/GSA/index.html)
**13.** [VennDiagram](https://cran.r-project.org/web/packages/VennDiagram/index.html)
**14.** [RColorBrewer](https://rdrr.io/cran/RColorBrewer/)
**15.** [ggVennDiagram](https://cran.r-project.org/web/packages/ggVennDiagram/index.html#:~:text=ggVennDiagram%3A%20A%20'ggplot2'%20Implement,geometry%20dataset%20and%20'ggplot2'.)


#### Bioconductor
**16.** [DEP](https://bioconductor.org/packages/release/bioc/html/DEP.html)
**17.** [vsn](https://www.bioconductor.org/packages/release/bioc/html/vsn.html)
**18.** [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
**19.** [BioNet](https://www.bioconductor.org/packages/release/bioc/html/BioNet.html)
**20.** [OmnipathR](https://bioconductor.org/packages/release/bioc/html/OmnipathR.html)
**21.** [mixOmics](https://bioconductor.org/packages/release/bioc/html/mixOmics.html)
**22.** [fgsea](http://bioconductor.org/packages/release/bioc/html/fgsea.html)


#### GitHub
**23.** [DEP](https://github.com/arnesmits/DEP)
**24.** [OmnipathR](https://github.com/saezlab/OmnipathR)
**25.** [NEMO](https://github.com/Shamir-Lab/NEMO)


To start a R-Console, you first need to make sure the R module is loaded

```bash
# check which modules are loaded
module list
# if the output contains e.g. R/4.0.5_deb10(default), you are fine...
# ...otherwise load the R module by typing "module load R"
# start R-Console
R
```

You can verify missing packages by running the following in the R-Console:

```python
package_list <- c("readr", "dplyr", "tidyr", "ggplot2", "ggrepel", "knitr", "tidyverse", "igraph", "ggpubr", "M2SMF", "SNFtool", "GSA", "VennDiagram", "RColorBrewer", "ggVennDiagram", "DEP", "vsn", "limma", "BioNet", "OmnipathR", "mixOmics", "fgsea", "NEMO")
print(setdiff(package_list, installed.packages()[, 1]))
```

You can install the missing packages by following the instructions below:

> Missing **CRAN** packages can be installed by simply running:

```python
install.packages("PackageName", dependencies = TRUE)
```

> Missing **GitHub** packages can be installed by simply running:

```python
require(devtools)
install_github('DeveloperName/PackageName')
```

> Missing **Bioconductor** packages can be installed by following the instructions provided on the Bioconductor web-pages dedicated to each of the packages or by running:

```python
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PackageName")
```
