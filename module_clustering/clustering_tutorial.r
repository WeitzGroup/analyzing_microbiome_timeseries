
## Library loading
library(reshape2); 
library(vegan); 
library(dplyr)
library(ade4); 
library(plotly)
library(compositions); 
library(pracma); 
library(DESeq2); 
library(fpc); 
library(tidyverse)
library(purrr)
library(cluster)
library(RColorBrewer)
library(ape)

# Migrate to location of repository and downloaded test data:
# setwd($PATH/analyzing_microbiome_timeseries/module1_clustering/)
load("TESTDATA_DIEL18S.RData",verbose=T)

# To address ecological questions, you may need to perform some quality 
# control on your OTU data. (1) Below, first global singletons are removed.
# This is unwanted noise in the form of OTUs which only have one sequence 
# in the entired dataset. 
# (2) Then, any OTUs assigned a Metazoan taxonomic identity are removed. This
# is because the methods used to generate this 18S data are not suitable
# to capture the metazoan community.
# (3) Finally, while we are mainly working with OTU.IDs throughout this tutorial
# we want to eventually add back in the taxonomic identities. So we make a taxonomy
# key dataframe ('tax_key')
#
# (1) Remove global singletons (where a single OTU appears with only 1 sequence in the entire dataset)
rowsum<-apply(test_data[2:20],1,sum) # sum numerics in OTU table by row
counts.no1 = test_data[ rowsum>1, ]  # remove rows that sum to 1 or less
dim(test_data)[1] - dim(counts.no1)[1] # report total number of singletons removed
#
# (2) Filter out OTUs which belong to the metazoa:
counts.filtered<-counts.no1[-(which(counts.no1$Level3 %in% "Metazoa")),] # Remove from counts
dim(counts.no1);dim(counts.filtered) # Test data should be left with 1984 total OTUs
#
# (3) Create taxonomy key:
names(counts.filtered)
seq_counts<-counts.filtered[1:20]
tax_key<-counts.filtered[c(1,21:31)]; head(tax_key[1:2,])
# Modify seq_counts data frame so row.names = OTU.ID and all values are numeric
row.names(seq_counts)<-seq_counts$OTU.ID
seq_counts$OTU.ID<-NULL; head(seq_counts[1:2,])

# Because these data are compositional, they are on a simplex. 
# The data being on a simplex means, if we know the abundances of all of the OTUs except for one, we still have the 
# same amount of information as if we knew all of the OTUs. 
# Evidence:
# Estimate covariance matrix for OTUs
covariance_matrix<-as.matrix(seq_counts)%*%t(seq_counts)
# Evaluate determinant of covariance matrix
cov_determinant<-det(covariance_matrix)
cov_determinant #0
# The determinant of the covariance matrix (what we just calculated) is equivalent to the product of the
# proportion of variance explained by every PCA axis. If the determinant is 0, that means there is an axis which
# explains 0 variance that we can't separate from the other axes. The data need to be transformed to be suitable
# for PCA. 

## Another approach to PCA ordination w/Compositional Data: Log-Ratio Transformations
# Log-ratio
log_rats<-data.frame(compositions::ilr(t(seq_counts)))
# These are complicated to interpret however because you move from D simplex to D-1 euclidean
# But using these we can see we are now in invertible covariance regime

## Checking the difference the log-ratio made on the data characteristics
new_covdet<-det(as.matrix(log_rats)%*%t(log_rats))
cov_determinant #Original Count Data
new_covdet # After doing a log ratio, you see that none of the axes are equal to 0

# AND change command so its the same.
lograt_pca<-prcomp(log_rats)
# Look at this
lograt_variances<-lograt_pca$sdev^2/sum(lograt_pca$sdev^2)
barplot(lograt_variances,
        main='Log-Ratio PCA Screeplot',
        xlab='PC Axis',
        ylab='% Variance',
       col=c(rep('black',1),rep('grey',18)),
       cex.names=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
legend('topright',fill=c('black','grey'),c('Should Present','??? Judgment Call'))

## We conclude a more faithful representation is to plot this ordination in 2D because the 3rd and 4th axes appear
## very similar, and we can't construct a 4D plot
pca_lograt_frame<-data.frame(lograt_pca$x,
                          time_of_day=gsub('^.*_','',rownames(lograt_pca$x)))
ggplot(pca_lograt_frame)+
    geom_point(aes(x=PC1,y=PC2,col=time_of_day))+
    ylab(paste0('PC2 ',round(lograt_variances[2]*100,2),'%'))+
    xlab(paste0('PC1 ',round(lograt_variances[1]*100,2),'%'))+
    scale_color_brewer(palette='Set1',name='Time of Day')+
    ggtitle('Log-Ratio PCA Ordination')+
    coord_fixed(ratio=lograt_variances[2]/lograt_variances[1])+
    theme_bw()

# So let's start by looking at PCoA. PCoA is doing a PCA on a distance matrix constructed from the data. 
# We then need a distance matrix. 
# Different distance metrics emphasize separate attributes/factors in microbial community comparison
# For instance, if we want to prioritize differences in presence/absence between samples
jac_dmat<-vegdist(t(seq_counts),method="jaccard") # p/a example
pcoa_jac<-ape::pcoa(jac_dmat) #perform PCoA
# Now we need to inspect how the PCoA turned out. We look at a Screeplot for this.
# The screeplot shows the amount of variance explained by each axis
samp_no<-dim(seq_counts)[2]
jac_variances<-pcoa_jac$values$Relative_eig
par(mar=c(5,6,4,1)+.1)
barplot(jac_variances,
        xlab='Principal Coordinate Axis',
        ylab='% Variance',
       col=c(rep('black',2),'darkgrey',rep('lightgrey',16)),
       cex.names=2,cex.axis=2,cex.lab=2,cex.main=1,
       names.arg=as.character(1:18))
legend('topright',fill=c('black','darkgrey','lightgrey'),c('Should Present','Judgment Call','Unnecessary to Present'),cex=2)
# How to read this plot: Before we plot the actual ordination, we need to decide which axes to present.
# We need to select as few axes as possible (so we can visualize) which capture large amounts of variance.
# In this example, we see the first axis captures most of the variance, so we definitely will show that one.
# Then the next two axes show less, but a similar amount, and the remaining all show way less in comparison. 
# Therefore, we can exclude everything after the 3rd axis. Because axes 2 and 3 capture similar amounts of variance,
# if we show one, we need to show the other one to be faithful to the data.


# Performing the log-ratio transformation makes the data all occupy a similar dynamic range, so we can use
# magnitude-sensitive distances like euclidean distance

euc_dmat<-dist(log_rats) 


## Conduct ordination w/distance matrix
pcoa_euc<-ape::pcoa(euc_dmat)
#euc_variances<-pcoa_euc$sdev^2/sum(pcoa_euc$sdev^2)
euc_variances<-pcoa_euc$values$Relative_eig
par(mar=c(5,6,4,1)+.1)
barplot(euc_variances,
        xlab='Principal Coordinate Axis',
        ylab='% Variance',
       col=c(rep('black',2),rep('darkgrey',2),rep('lightgrey',17)),
       cex.names=2,cex.axis=2,cex.lab=2,cex.main=1,
       names.arg=as.character(1:18))
legend('topright',fill=c('black','darkgrey','lightgrey'),c('Should Present','Judgment Call','Unnecessary'),cex=2)
#
# By trying 2 different metrics, we see a differences in ordination output.
# For Jaccard distance, we see we need 3 axes to present the ordination. Using euclidean distance,
# because of the complex covariance structure of the relative compositions, the data do not easily
# ordinate into 2 or 3 dimensions.


# How dissimilar are my samples with respect to diversity?
## P/A is recommended when asking questions about the ENTIRE community
# Does the active microbial community change depending on time of day?
# To emphasize changes in presence/absence, we elect to use Jaccard distance for ordination.
# Our ordination indicates the data are most suitably summarized
# with 3 dimensions. 
#pcoa_jac_frame<-data.frame(pcoa_jac$x,
#                           time_of_day=gsub('^.*_','',rownames(pcoa_jac$x)))
pcoa_jac_frame<-data.frame(pcoa_jac$vectors,time_of_day=gsub('^.*_','',rownames(pcoa_jac$vectors)))
eigenvalues<-round(jac_variances,4)*100
plot_ly(pcoa_jac_frame,type='scatter3d',mode='markers',
        x=~Axis.2,y=~Axis.3,z=~Axis.1,colors=~brewer.pal(6,'Set1'),color=~time_of_day)%>%
  layout(font=list(size=18),
         #title='PCoA Jaccard Distance',
         scene=list(xaxis=list(title=paste0('Co 2 ',eigenvalues[2],'%'),
                    showticklabels=FALSE,zerolinecolor='black'),
         yaxis=list(title=paste0('Co 3 ',eigenvalues[3],'%'),
                    showticklabels=FALSE,zerolinecolor='black'),
         zaxis=list(title=paste0('Co 1 ',eigenvalues[1],'%'),
                    showticklabels=FALSE,zerolinecolor='black'),
        aspectratio = list(x=3*eigenvalues[2]/eigenvalues[1])), y=3*eigenvalues[3]/eigenvalues[1], z=3)

## To compare with the euclidean distance ordination representing differences in relative composition
#pcoa_euc_frame<-data.frame(pcoa_euc$x,
#                          time_of_day=gsub('^.*_','',rownames(pcoa_euc$x)))
pcoa_euc_frame<-data.frame(pcoa_euc$vectors,time_of_day=gsub('^.*_','',rownames(pcoa_euc$vectors)))
euc_eigenvalues<-round(euc_variances,4)*100
plot_ly(pcoa_euc_frame,type='scatter3d',mode='markers',
        x=~Axis.3,y=~Axis.2,z=~Axis.1,colors=~brewer.pal(6,'Set1'),color=~time_of_day)%>%
  layout(font=list(size=18),
         #title='PCoA Euclidean Distance',
         scene=list(xaxis=list(title=paste0('Co 3 ',euc_eigenvalues[3],'%'),
                    showticklabels=FALSE,zerolinecolor='black'),
         yaxis=list(title=paste0('Co 2 ',euc_eigenvalues[2],'%'),
                    showticklabels=FALSE,zerolinecolor='black'),
         zaxis=list(title=paste0('Co 1 ',euc_eigenvalues[1],'%'),
                    showticklabels=FALSE,zerolinecolor='black'),
                   aspectratio = list(x=3*euc_eigenvalues[3]/euc_eigenvalues[1], 
                                      y=3*euc_eigenvalues[2]/euc_eigenvalues[1], 
                                      z=3)))
## Further note: If your ordination has data which align in a 'T' or '+' shape perpendicular to the axes
## this is often diagnostic of covariance attributed to the higher dimensions which are not plotted

# To corroborate the 3-D plot, use a simple average 
## hierarchical clustering and plot the dendrogram.
## Some of the same pattern emerge, but it is not as 
## clear as the 3-D representation.
cluster_ex<-hclust(vegdist(t(seq_counts),method='jaccard'),method="complete") #Using Jaccard distance so apples to apples
plot(cluster_ex,main='Jaccard Hierarchical Agglomerative Clustering',xlab='',sub='')

# Although our statistical ordinations appear to require at least 3 dimensions to communicate the data.
# However, we don't always have the budget for a 3D plot. So we may want to impose the condition on an ordination
# technique that the answer MUST go in 2D. We turn to NMDS here. 
#
set.seed(071510) # setting random seed to assure NMDS comes out the same every time
## So we can compare the relative composition based distance metric to the presence/absence based distance metric
euc_nmds<-metaMDS(euc_dmat,k=2,autotransform=FALSE)
jac_nmds<-metaMDS(jac_dmat,k=2,autotransform=FALSE)

# Take a look at stress - overall this value is not extremely informative, but know that the
# closer stress is to 1 the less representative of your actual data the NMDS is
euc_nmds$stress #Stress ~0.075
jac_nmds$stress #Stress ~0.073 So the Jaccard NMDS is a ***slightly*** more parsimonious ordination.

# Additionally, the axes for NMDS are totally arbitrary, so axis scaling does not matter
# and data can be rotated/reflected about axes and the NMDS is still the same
euc_frame<-data.frame(euc_nmds$points,
                         time_of_day=gsub('^.*_','',rownames(log_rats)))
jac_frame<-data.frame(jac_nmds$points,
                     time_of_day=gsub('^.*_','',rownames(log_rats)))
## Plotting euclidean distance NMDS
ggplot(euc_frame,aes(x=MDS1,y=MDS2,col=time_of_day))+
    geom_point(size=2)+
    scale_color_brewer(palette='Set1',name='Time of Day')+
    theme_bw()+ggtitle('Euclidean Distance NMDS')
## Plotting Jaccard distance NMDS
ggplot(jac_frame,aes(x=MDS1,y=MDS2,col=time_of_day))+
    geom_point(size=2)+
        scale_color_brewer(palette='Set1',name='Time of Day')+
    theme_bw()+ggtitle('Jaccard Distance NMDS')
## For paper figures
euc_fig<-ggplot(euc_frame,aes(x=MDS1,y=MDS2,col=time_of_day))+
    geom_point(size=10)+
    scale_color_brewer(palette='Set1',name='Time of Day')+
    theme_bw()+theme(legend.position='none',text=element_text(size=24))+
xlab('NMDS 1')+ylab('NMDS 2')
jac_fig<-ggplot(jac_frame,aes(x=MDS1,y=MDS2,col=time_of_day))+
    geom_point(size=10)+
    scale_color_brewer(palette='Set1',name='Time of Day')+
    theme_bw()+theme(legend.position='none',text=element_text(size=24))+
xlab('NMDS1')+ylab('NMDS2')
euc_fig
jac_fig

## Before we can cluster temporal dynamics, we must make sure species are intercomparable. 
## We need to perform normalization steps to ensure intercomparability between species.
## Put another way, the count data are heteroskedastic -- 
## the mean is correlated with the variance.
within_seq_means<-apply(t(seq_counts),2,mean)
within_seq_vars<-apply(t(seq_counts),2,var)
plot(within_seq_means,within_seq_vars,log='xy',
    main='Heteroskedastic Data',
    xlab='Mean # Counts',
    ylab='Var # Counts')
# Variance (in the amount of times an OTU shows up) is correlated to average count. 
# This attribute (common to count data) violates the assumption of constant variance in errors
# which is necessary for standard linear modeling. To use z-scores for comparisons and to detrend our timeseries 
#(a kind of linear model), we have to deal with this.

## We stabilize the variances so that they are no longer correlated to the species abundances.
transformed_counts<-DESeq2::varianceStabilizingTransformation(as.matrix(seq_counts))

## Check to see how the transformation did
within_trans_means<-apply(transformed_counts,1,mean)
within_trans_vars<-apply(transformed_counts,1,var)
# Plot:
plot(within_trans_means,within_trans_vars)
# See how now most of the variances are the same between otus? 
# This means the data are more intercomparable using a z-score transformation + detrending

# Since the goal is to identify repeating day/night patterns, we need to remove trends we are not considering,
# i.e., linear trends
#?pracma::detrend() for more information
transformed_detrended<-apply(transformed_counts,1,pracma::detrend)
# Now we rescale so numbers are again intercomparable and now overdispersion has been
# dealt with
trans_dt_scaled<-apply(transformed_detrended,2,scale)
## We lost rownames so tack those back on
rownames(trans_dt_scaled)<-colnames(transformed_counts)
#
## Comparing the distribution of the data along every step of the pipeline -- notice how by the time we get to the 
## z-scores we're looking more like a standard normal
hist(transformed_counts,
    main='VST Sequence Count Observations',
    xlab='# Total Observations',
    ylab='Frequency')
hist(transformed_detrended,
    main='VST+Detrended Data',
    xlab='Total Observations Accounting for Linear Trends',
    ylab='Frequency')
hist(trans_dt_scaled,
    main='VST+Detrended+Scaled Data (Z-scores)',
    xlab='Expression Level Relative to per-OTU Avg',
    ylab='Frequency')

# Now we're ready to generate a distance matrix for clustering
# Recall what we previously learned about different distance metrics.
# We use euclidean distance below:
temporal_dmat<-dist(t(trans_dt_scaled))
n_clusts<-2:20 # This sets the range of # clusters we divide the data into
# First, the classic hierarchical agglomerative clustering
hc_full_cluster<-hclust(temporal_dmat)
hc_clusts<-lapply(n_clusts,function(x) cutree(hc_full_cluster,x))
## Now we try k-medoids, also a very popular clustering method
## This particular implementation is not optimized for speed, so this may take a minute or so
kmed_clusts<-lapply(n_clusts, function(x) cluster::pam(temporal_dmat,k=x))

## Comparing cluster outcomes
hc_stats<-lapply(hc_clusts,function(x) fpc::cluster.stats(temporal_dmat,
                                                          clustering=x))
kmed_stats<-lapply(kmed_clusts, function(x) fpc::cluster.stats(temporal_dmat,
                                                             clustering=x$clustering))

## We write a helper function to be less redundant
ripping_stats<-function(list,func){
  ## Essentially all this function does is implements a function (func) on a list
  ## and coerces the output to a column vector (this will be handy when we want to make a data frame)
  output<-do.call(rbind,lapply(list,func))
  return(output)
}
func_list<-rep(list(function(x) x$cluster.number,
                function(x) x$within.cluster.ss,
                function(x) x$avg.silwidth,
                function(x) x$ch),2)
stats_list<-rep(list(hc_stats,kmed_stats),each=4)
collected_stats<-purrr::map2(stats_list,func_list,ripping_stats)
nclusts<-rep(n_clusts,length(collected_stats))
method<-rep(c('hc_agg','kmed'),each=length(n_clusts)*length(collected_stats)/2)
ind_name<-rep(c('n','ss','sil','ch'),each=length(n_clusts)) %>%
  rep(2)
                    
## Collecting outputs for plotting                    
index_frame<-data.frame(index=do.call(rbind,collected_stats),
                        nc=nclusts,
                        Method=method,
                        ind=ind_name)
index_frame %>%
  filter(ind=='ss') %>%
  ggplot(aes(x=nc,y=index,col=Method)) +
  geom_point() +
  geom_line(aes(group=Method)) +
  ylab('Within Cluster Sum Square Error') +
  xlab('Number Clusters')+
  theme_bw()
## Interpretation of plots:                    
# As we add more clusters (yaxis), we want the sum square error to decrease
# As we allow for more clusters, we are able to describe more of the variability in the data
# so the sum square error decreases.
## trade-off - it will always decrease when you add clusters, but we want to 
## add as fewer clusters as possible, hence in this example kmed > hc
# for any given amount of clusters, kmed has less error. So we will continue using kmed approach.   

## Plotting out the calinski-harabasz indices for each clustering
## We want this number to be high
index_frame %>%
  filter(ind=='ch') %>%
  ggplot(aes(x=nc,y=index,col=Method)) +
  geom_point() +
  geom_line(aes(group=Method)) +
  ylab('C-H Stat') +
  xlab('Number Clusters')+
theme_bw()

## Silhouette profile: Looking at species variability w/in our 8 clusters
## Extract silhouette profile for the clustering we select
silhouette_profile_kmed8<-cluster::silhouette(kmed_clusts[[7]])
## Reformat for easier plotting
silhouette_frame<-data.frame(cluster=silhouette_profile_kmed8[,1],
                             neighbor=silhouette_profile_kmed8[,2],
                             dist=silhouette_profile_kmed8[,3],
                             otu_id=rownames(silhouette_profile_kmed8))

## Sorting by cluster and ranking by silhouette width
new_silframe<-silhouette_frame %>%
  arrange(dist) %>%
  group_by(cluster) %>%
  mutate(idno=1:n(),
         tot_num=n())

## Plotting silhouette profile of OTUs for each cluster
ggplot(new_silframe,aes(x=idno,y=dist))+
  geom_bar(stat='identity') +
  coord_flip() +
  ggtitle('Silhouette Profiles for Individual Clusters') +
  facet_wrap(~cluster)

## In order to interpret the silhouette profiles, we can compare the width (x-axis) and height (y-axis) of each
## profile. For instance, comparing cluster #3 vs #8: Cluster 3 is much taller (has more species in it) compared to 
# cluster 8, but also is much broader (there are more species with intermediate values on the x-axis). This means
# on the whole cluster 3 is less coherent than cluster 8. 
## Also keep an eye out for species w/negative silhouette distances -- this means that they have a nearby neighbor
## who belongs in a different cluster and therefore may be only marginally similar to their cluster on the whole

# What is the taxonomic composition of each of the 8 clusters?
## To make ecological sense of the temporal clusters,
## we can pull out a representative species and use their signal over time
## to summarize the temporal trends across the whole cluster
medoid_otus<-kmed_clusts[[7]]$medoids
## Extracting those time series
medoid_dynamics<-trans_dt_scaled[,medoid_otus]

# Using the output from:
head(silhouette_frame)
# We can re-build what the taxonomic composition of each 
## of these clusters is.
range(silhouette_frame$cluster) # n=8 clusters with test data

# Join cluster information with taxonomic information
colnames(tax_key)[1]<-"otu_id"
tax_bycluster<-left_join(silhouette_frame, tax_key, by="otu_id") # Recall the taxonomy key we made from above?
head(tax_bycluster[1:2,])

# Summarize taxa with cluster information
tax_bycluster_sum<-tax_bycluster %>%
    group_by(cluster, Taxa) %>%
    summarise(richness =  n()) %>% #Get richness information
    group_by(Taxa) %>%
    mutate(n_per_tax=sum(richness)) %>% #Figure out total number of OTUs assigned to each taxon
    as.data.frame
tax_bycluster_sum$Taxa[which(tax_bycluster_sum$Taxa=='Opisthokonts')]<-'Opisthokont' # Cleaning this up
head(tax_bycluster_sum)

# Plot cluster-to-cluster taxonomic differences
## Values are equal to the total number of OTUs 
## belonging to each taxonomic group
#
tax_color<-c("#67000d","#e31a1c","#dd3497","#fcbba1","#fed976","#fc8d59","#a63603","#addd8e","#7f2704","#238b45","#a1d99b","#081d58","#1f78b4","#a6cee3","#8c6bb1","#9e9ac8","#984ea3","#081d58","#662506","#ffffff","#969696","#525252","#000000")
#
ggplot(tax_bycluster_sum, aes(x=cluster, y=richness, fill=Taxa))+
    geom_bar(color="black", stat="identity")+
    coord_flip()+
    scale_fill_brewer(palette='Set3')+
    theme_minimal()+
    scale_x_continuous(breaks=1:8,labels=as.character(1:8),name='Cluster #')+
    scale_y_continuous(name='# Species')+
theme(legend.position='bottom',text=element_text(size=18))+ ## Make the cluster #s happen
guides(fill=guide_legend(nrow=3,byrow=TRUE))

## Demonstrating medoid dynamics for each cluster
medoid_dyn_long<-as.data.frame(t(medoid_dynamics)) %>%
mutate(otu_id=colnames(medoid_dynamics),
      clust_num=1:8) %>%
gather(timepoint,z_score,'1_6PM':'19_6PM') %>%
mutate(time_numeric=as.numeric(gsub('_.*$','',timepoint)),
      plotting_timepoint=paste('Day',ceiling(time_numeric/6),gsub('^.*_','',timepoint))) 
#
## We're plotting the y axis as z-scores so that transcripts with different baseline levels are 
# more readily intercomparable. Read a z-score y-axis this like this:
# 0 is the average transcript level
# negative numbers are below-average numbers of transcripts
# positive numbers are above-average levels of transcripts
# 1 unit is 1 standard deviation away from the mean (So a value of 1 corresponds to mean+1sd, -1 is mean-1sd)

ggplot()+
geom_point(size=4,
           shape=21, 
           color="white", 
           aes(fill=factor(clust_num),
              x=time_numeric,
              y=z_score),
          data=medoid_dyn_long)+
geom_line(data=medoid_dyn_long,
          aes(x=time_numeric,
              y=z_score,
              group=otu_id,
              col=factor(clust_num)),size=1.25,linetype=5)+
facet_wrap(~clust_num,ncol=2)+
scale_color_brewer(palette='Set2',name='Cluster #')+
scale_fill_brewer(palette='Set2',guide=FALSE)+
scale_x_continuous(breaks=seq(1,19,by=6),labels=unique(medoid_dyn_long$plotting_timepoint)[c(1,7,13,19)])+
scale_y_continuous(name='Z-score Expression',limits=c(-2,4.5))+
theme(axis.text.x=element_text(angle=90,hjust=0.5,size=18),
      panel.background=element_rect(fill='white'),
     panel.grid.major=element_line(color='gray'),
     axis.title.x=element_blank(),
     text=element_text(size=18))+
geom_rect(data=data.frame(ymin=rep(-2,4),
          ymax=rep(4.5,4),
          xmin=seq(1,19,by=6),
         xmax=c(seq(4,19,by=6),19)),
          aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax),
          col='gray',
          alpha=0.25)
