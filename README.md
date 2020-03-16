## **A Primer for Microbiome Time-Series Analysis**

Code tutorials that accompany Coenen et al. _In Review_. 
Citation: Ashley R Coenen, Sarah K Hu, Elaine Luo, Daniel Muratore, and Joshua S Weitz. *A Primer for Microbiome Time-Series Analysis* _In Review_.  

### Summary  
Time-series can provide critical insights into the structure and function of microbial communities. The analysis of temporal data warrants statistical considerations, distinct from comparative microbiome studies, to address ecological questions. This primer identifies unique challenges and approaches for analyzing microbiome time-series. In doing so, we focus on (1) identifying compositionally similar samples, (2) inferring putative interactions among populations, and (3) detecting periodic signals. 

These tutorials are meant to serve as a primer that includes the basics of how to approach statistical ordination, clustering, parameter estimation, and nonparametric periodic signals in high-resolution temporal data.  

#### **Tutorials**
* Introduction (**MATLAB, Octave**) - demonstrates how autocorrelation within timeseries leads to spurious correlations between independent variables
* Clustering (**R**) - uses a test dataset to explore normalization of compositional data, diagnostic tools to ensure ordination approaches are appropriate, and ordination methods including NMDS and PCoA
* Periodicity (**R**) - with a test dataset, demonstrates normalization of compositional data and how to pull out periodic signals 
* Inference (**MATLAB**) - with simulated data, demonstrates how to perform model-based network inference to estimate the infection network in a microbe-virus community
