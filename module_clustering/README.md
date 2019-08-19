# Introduction for running module on Statistical Ordination and Clustering
## Section 1: Getting started
_Setting up your working environment_ 

There are two ways to run the provided code, first is using **jupyter notebooks with R** or open the provided code in **RStudio**. An example of how to run RStudio and jupyter notebooks using conda environments can be found [here](https://alexanderlabwhoi.github.io/post/anaconda-r-sarah/).  

### jupyter notebooks

Activate an R conda environment, migrate to this repo, and start jupyter notebooks to run ```clustering_tutorial.ipynb```.

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

To load the test data set:
```
load("TESTDATA_DIEL18S.RData",verbose=T)
```
This data is an OTU table with taxonomy, where each row is an OTU and columns are samples. The last columns in the data list the taxonomic identity.  

**Tutorial workflow:**
1. Import tag sequence dataset
   - Remove global singletons
   - Filter out unwanted OTUs
   - Generate a separate taxonomy dataframe so we can work with all numeric data
2. Estimate covariance matrix for OTUs 
3. Perform log-ratio transformation
   - Evaluate how the transformation altered the covariance
   - Perform principal component analysis on log-ratio transformed data
   - Plot resulting PC axes as a bar plot
   - Plot PCA in 2 dimentions
4. Calculate a Jaccard distance matrix from the count data
   - Perform a PCoA with the resulting Jaccard dist matrix
   - Plot PC axes as a bar plot
5. Calculate Euclidean distance from log-ratio transformed data
   - Perform PCoA
   - Plot PC axes as a bar plot
6. Plot 3 dimensional PCoAs using plot_ly
   - Plot 3-D Jaccard PCoA
   - Plot 3-D Euclidean PCoA
   - Compare results with a Jaccard Hierarchical Agglomerative cluster plot
7. Calculate Nonmetric multidimensional scaling
   - Calculate with both Euclidean and Jaccard transformed data
   - Compare stresses of each
   - Plot 2 dimensions MDS plot with Euclidean and Jaccard transformed data
8. Address 'Are there taxa within samples that co-occur?'
   - Demontrate data are heteroskedastic
   - Perform Variance Stabilizing Transformation using DESeq
   - Evaluate how this altered the mean and variance (plot)
9. De-trend transformed count data
   - Compare the distribution of transformed count data before and after detrending and scaling
10. Generate a distance matrix from transformed data
   - Use hierarchical agglomerative clustering
   - Use k-medoids clustering
   - Compare clustering approaches
   - Plot Calinski-Harabasz indices from each cluster
11. Generate a Silhouette profile to evaluate clustering approach
   - Dissect clustered data by taxonomic identity
   - Demonstrate medoid dynamics for each cluster

***

If you have questions please feel free to reach the authors at the contact provided in the main README
