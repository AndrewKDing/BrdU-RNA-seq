library(readxl)
library(tidyverse)
library(ggplot)

# CHECKPOINTS
# Calculate MDS and plot in ggplot
      # CHECKPOINT: Calculate 
# PRIORITY PLOT: colored by date of library preparation
# PRIORITY PLOT: colored by date of RNA extraction
# PRIORITY PLOT: colored by date of sorting
# PLOT: BrdU+, colored by non-infected, 21dpi, and necropsy
# PLOT: BrdU-, colored by non-infected, 21dpi, and necropsy
# PLOT: BrdU+, colored by %BrdU
# PLOT: 21dpi only, non-infected, SIVE, SIVnoE
# plot: necropsy only, non-infected, SIVE, SIVnoE


# Principal components analysis -------------------------------------------

#start with prcomp, you need to transpose your data first, because prcomp takes
#samples as rows

# need to multiply lib.size by norm.factors, or just use getCounts?
# I need to calculate the normalized counts.. plotMDS uses logCPM values but
# that's not how dream is normalizing

#function one, returns principal components of dgelist
pcaobj <- prcomp(t(dge[[1]]))
pca_norm <- cpm(dge, log=TRUE, prior.count=3)
pca_norm_obj <- prcomp(t(pca_norm))

#pulling out the principal components
pca_matrix <- as.data.frame(pca_norm_obj$x)

#the "x" variable holds the principal components
#you need to merge in the sample metadata here. merge by the 0th column (sample id)


barplot(round(pca_norm_obj$sdev^2/sum(pca_norm_obj$sdev^2)*100,2),las=2,
         names.arg=colnames(pca_norm_obj$x),ylab="% Variance explained",
         xlab="PCA principal components")
#Alright, after normalization it's still different. First principal component is
#50% of variance. Is something wrong with plotMDS?

#this is, oddly, giving me a different result where the first principal component explains essentially all of the variance. what gives?
plotExpressionPCA(dge)
#https://support.bioconductor.org/p/103747/#103769 source on normalization

#going fro 500 to 20,000 genes matters a LOT. it increases principle component 1 variance from like 10% to 42%. I wonder if this means something about 
# gene signatures instead of specific differentially expressed genes?
# Using more genes looks correct, non-infected samples appear to be clustering with one another



#Using dist followed by cmdcale to plot MDS

#Note to self, data should be subsetted PRIOR to multi-dimensional scaling

distances <- cmdscale(dge[[1]])
test <- as.matrix(count_data[,-(1:2)])
rownames(test) <- paste0(count_data[,1], "=",count_data[,2])

plotMDS(dge)


#should we remove batch effects?
#let's make a number of PCA plots
#we should do this in ggplot
#pl