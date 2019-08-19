
library(compositions)
library(pracma)
library(rain)
library(dplyr)
library(reshape2)
library(ggplot2)

# Migrate working directory to the location of repository and downloaded test data
# setwd() assuming data is in folder with script
load("TESTDATA_DIEL18S.RData",verbose=T) # Import test data set
head(test_data)

# Repeat from Module 1 Clustering tutorial (QC initial OTU table)
# (1) Remove global singletons:
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

# Centered log-ratio transformation. 
# This is necessary due to the compositional nature of OTU results 
# ?clr()
df_clr<-as.data.frame(clr(seq_counts))

# Detrend
# ?detrend()
df_detr<-detrend(t(df_clr))

# RAIN analysis
# ?rain()
# Set up RAIN parameters. 
# This analysis will treat each time point independently (rather than each 'day' as a replicate)
t<-(1:19)
ft <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

# Run RAIN
output_RAIN<-rain(as.matrix(df_detr), period=24, measure.sequence=ft, deltat=4, method="independent", na.rm = TRUE)

# Report results
head(output_RAIN[1:3,]) # Output data frame

hist(output_RAIN$pVal) # distribution of p-value results

sig<-subset(output_RAIN, pVal < 0.05); dim(sig)[1] # Total number of OTUs with significant rhymicity (p<0.05)

sig$OTU.ID<-row.names(sig); list_of_sig_OTUs<-unique(sig$OTU.ID) # list of OTUs with significant rhythmicity

# Assessing significance
# Bonferroni correction
# This is the LEAST forgiving form of multiple testing correction (least power and least type I error)
bonf_sig<-subset(output_RAIN, pVal <= (0.05/nrow(output_RAIN)))

# Benjamini-Hochberg Correction
#?p.adjust()
fdr_ps<-p.adjust(output_RAIN$pVal,method='BH')
bh_sig<-output_RAIN[which(fdr_ps<=0.05),]

# Adaptive Benjamini-Hochberg Correction (see Benjamini and Hochberg 2001)
# Step 1: Sort p-values in ascending order
abh_ps<-output_RAIN$pVal[order(output_RAIN$pVal,decreasing=FALSE)]

# Step 2: Find p-values smaller than those expected at 5% FDR.
# Rejecting all null hypotheses which have p-values less than their corresponding abh_metric will
# yield the same results as p.adjust(method='BH') at the 0.05 significance level.
q<-0.05 #Significance level
m<-length(abh_ps) #Number of hypotheses tested
abh_metric<-q*(1:m)/m

#For a visualization of what this procedure does -- we consider all red points that fall below the black
# line significant
plot(1:m,abh_metric,type='l',
  xlab='p-value rank',
  ylab='p-value',lwd=2,
  main='Benjamini-Hochberg FDR Control')
points(1:m,abh_ps,col='red')
legend(650,0.05,fill=c('red','black'),legend=c('p-value','FDR controlled p threshold'))

# Step 3: Estimate the number of false null hypotheses according to the procedure outlined
# in Benjamini and Hochberg 2001
s=0
s_i=0
i=1
while(s>=s_i){
  s_i=s
  s=(1-abh_ps[i])/(1+m-i)
  i=i+1
}
m_0<-ceiling(1/s)+1

# Step 4: Find p-values smaller than those expected at 5% FDR for adjusted number of false hypotheses
abh_updated<-(1:length(abh_ps))*0.05/m_0
# For a similar visualization
plot(1:m,abh_updated,type='l',
     xlab='p-value rank',
     ylab='p-value',lwd=2,
     main='Adaptive Benjamini-Hochberg FDR Control')
points(1:m,abh_ps,col='red')
legend(650,0.05,fill=c('red','black'),legend=c('p-value','FDR controlled p threshold'))

# Step 5: Determine which tests resulted in a p-value significant at the 5% level considering FDR
abh_sigs<-output_RAIN[which(output_RAIN$pVal %in% abh_ps[which(abh_ps<=abh_updated)]),]

# Compile list of significantly perodic OTUs with taxonomic identities (use tax_key generated above)
abh_sigs$OTU.ID<-row.names(abh_sigs)
sig_bytax<-left_join(abh_sigs, tax_key, by="OTU.ID")
dim(abh_sigs); dim(sig_bytax)
# Report list of OTUs found to have significant diel rhythmicity
write.table(sig_bytax, file="SigDiel_OTUs.txt", row.names = FALSE, quote = FALSE)

# Join CLR normalized data with list of significantly periodic OTUs
df_clr2<-df_clr
df_clr2$OTU.ID<-row.names(df_clr)
df_clr3<-left_join(sig_bytax, df_clr2, by="OTU.ID")
unique(df_clr3$Taxa) # What is the overal taxonomic distribution?
head(df_clr3)

# Plot signal of OTUs with significant periodicity
df_clr4<-melt(data.frame(df_clr3$OTU.ID, df_clr3$Taxa, df_clr3$taxonomy, df_clr3[17:35])) #Extract taxa names and CLR normalied values
df_clr4$Time_of_day<-factor(df_clr4$variable, 
                            levels = c("X1_6PM","X2_10PM","X3_2AM","X4_6AM","X5_10AM","X6_2PM","X7_6PM","X8_10PM","X9_2AM","X10_6AM","X11_10AM","X12_2PM","X13_6PM","X14_10PM","X15_2AM","X16_6AM","X17_10AM","X18_2PM","X19_6PM"), 
                            labels = c("Day 1 6PM","Day 1 10PM","Day 2 2AM","Day 2 6AM","Day 2 10AM","Day 2 2PM","Day 2 6PM","Day 2 10PM","Day 3 2AM","Day 3 6AM","Day 3 10AM","Day 3 2PM","Day 3 6PM","Day 3 10PM","Day 4 2AM","Day 4 6AM","Day 4 10AM","Day 4 2PM","Day 4 6PM"))
colnames(df_clr4)[1:3]<-c("OTU.ID","Taxa","Full taxonomy")
head(df_clr4)
#
# Plot CLR normalized values for those OTUs with significant diel rhythmicity"
diel_plot<-ggplot(df_clr4, aes(x=Time_of_day, y=value, group= OTU.ID))+
    geom_point(size=4,shape=21, color="white", aes(fill=Taxa, group = OTU.ID))+
    geom_line(size=1, linetype = 5, aes(color=Taxa))+
    facet_wrap(OTU.ID ~ `Full taxonomy`,labeller = label_wrap_gen(40))+
    theme_minimal()+
    labs(title="Taxa with significant diel rhythmicity", x = "Time of Day", y = "CLR normalized value")+
    theme(legend.position="none",axis.text.x = element_text(angle=45,hjust=1, vjust=1), 
          strip.text = element_text(hjust=0, size=8))+
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=4,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+
    geom_rect(data=NULL,aes(xmin=7,xmax=10,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+
    geom_rect(data=NULL,aes(xmin=13,xmax=16,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+
    geom_rect(data=NULL,aes(xmin=19,xmax=Inf,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)


diel_plot<-ggplot(df_clr4, aes(x=Time_of_day, y=value, group= OTU.ID))+
    geom_point(size=4,shape=21, color="white", aes(fill=Taxa, group = OTU.ID))+
    geom_line(size=1, linetype = 5, aes(color=Taxa))+
    facet_wrap(~`Full taxonomy`,labeller = label_wrap_gen(40))+
    theme_minimal()+
    labs(x = "Time of Day", y = "CLR normalized value")+
    theme(legend.position="none",axis.text.x = element_text(angle=70,hjust=1, vjust=1,size=12), 
          strip.text = element_text(hjust=0, size=12),axis.title=element_text(size=14))+
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=4,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+
    geom_rect(data=NULL,aes(xmin=7,xmax=10,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+
    geom_rect(data=NULL,aes(xmin=13,xmax=16,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+
    geom_rect(data=NULL,aes(xmin=19,xmax=Inf,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)
#
# Plot a subset
tax_of_interest<-c("Stramenopiles", "Haptophytes")
diel_plot %+% subset(df_clr4, Taxa %in% tax_of_interest)+
    scale_fill_manual(values=c("#de77ae", "#7fbc41"))+
    scale_color_manual(values=c("#de77ae", "#7fbc41"))

# Plot all
diel_plot


