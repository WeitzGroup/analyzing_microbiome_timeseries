# Introduction for running module on Statistical Ordination and Clustering

## Section 1: Getting started
_Setting up your working environment_ 

There are two ways to run the provided code, first is using **jupyter notebooks with R** or open the provided code in **RStudio**. An example of how to run RStudio and jupyter notebooks using conda environments can be found [here](https://alexanderlabwhoi.github.io/post/anaconda-r-sarah/).  

Another reason we recommend running this tutorial in a conda environment, is to streamline the R package installations.

### Here are the R packages that need to be installed ahead of time:

> reshape2 , vegan , dplyr, ade4 , plotly, compositions , pracma , DESeq2 , fpc , tidyverse, purrr, cluster, RColorBrewer, ape

To install using conda most packages are already available in the r-essentials package
```
conda install -c r r-essentials
```
All other packages can be installed by searching for conda install options, some of which may require a bioconductor install:
```
conda install -c bioconda bioconductor-deseq
```

To install in R command line or RStudio, insert package name in quotes to install:
```
install.packages("reshape2")
```
For those that require bioconductor, see instructions [online](https://www.bioconductor.org/install/):
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```
* Here is how to load a conda environment/how to install conda
* Here is how to install an R package via conda-install bioconductor [link_to_bioconductor_bioconda_website] 
* Here are the packages you need

# Section 2: Running the code
* Here is how to use jupyter notebooks


If you have questions please feel free to reach the authors at the contact provided in the main README
# Other packages can be installed using bioconductor
# Other packages can be installed using bioconductor
# Other packages can be installed using bioconductor
