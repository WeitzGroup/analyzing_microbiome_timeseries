# Introduction for running module on Statistical Ordination and Clustering

## Section 1: Getting started
_Setting up your working environment_ 

There are two ways to run the provided code, first is using **jupyter notebooks with R** or open the provided code in **RStudio**. An example of how to run RStudio and jupyter notebooks using conda environments can be found [here](https://alexanderlabwhoi.github.io/post/anaconda-r-sarah/).  

### jupyter notebooks

Activate an R conda environment, migrate to this repo, and start jupyter notebooks to run ```clustering_tutorial.ipynb```. To run$

```
# Example
conda activate r_3.5.1

cd ../../../analyzing_microbiome_timeseries/module_clustering/

jupyter notebook
# Open clustering_tutorial.ipynb
```

### R command line or RStudio
Open terminal and active R or start RStudio to run ```clustering_tutorial.r```.

## Section 2: R package installation

Here are the R packages that need to be installed ahead of time:
> reshape2 , vegan , dplyr, ade4 , plotly, compositions , pracma , DESeq2 , fpc , tidyverse, purrr, cluster, RColorBrewer, ape

To install using conda most packages are already available in the r-essentials package. *Using conda environments is also recommended, as R package installations may be easier than installing directly in RStudio or on the command line.
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


## Section 3: Running the code

Once your R environment is set up you can begin running the tutorial. The tutorial includes comments that parallel the primer. The included test data are from [Hu et al. 2018](https://www.frontiersin.org/articles/10.3389/fmars.2018.00351/full) additional information can be found [here](https://github.com/shu251/18Sdiversity_diel).



If you have questions please feel free to reach the authors at the contact provided in the main README
