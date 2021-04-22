
# *methcp*: An R package for Differentially Methylated Region Detection

*MethCP* is an R package for analyzing whole-genome bisulfite sequencing (WGBS) data. It implements *MethCP* algorithm introduced in [(Gong and Purdom, 2018)](https://www.biorxiv.org/content/early/2018/02/13/265116).

*MethCP* detects differentially methylated region for WGBS data based on change point detection, which naturally segments the genome and provide region-level differential analysis. It uses as input the results of a per-nucleotide test-statistic, like [DSS](https://www.ncbi.nlm.nih.gov/pubmed/24561809) or [methylKit](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r87), and uses this input to segment the genome into regions. Then, several meta-analysis strategies are provided to identify which of those regions are DMR. 

To install the package: 

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MethCP")
```

Alternatively, to install the current version on Github:

```r
library(devtools)
install_github("boyinggong/methcp", dependencies=TRUE,
  repos=BiocInstaller::biocinstallRepos())
```

A tutorial can be found [here](http://www.bioconductor.org/packages/release/bioc/vignettes/MethCP/inst/doc/methcp.html).

Contact: Boying Gong([boyinggong@berkeley.edu](mailto:boyinggong@berkeley.edu)).
