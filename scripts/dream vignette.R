library("variancePartition")
library("edgeR")
library("BiocParallel")
library("tidyverse")

data(varPartDEdata)

# filter genes by number of counts
isexpr <- rowSums(cpm(countMatrix) > 0.1) >= 5

# Standard usage of limma/voom
dge <- DGEList(countMatrix[isexpr, ])
dge <- calcNormFactors(dge)

# make this vignette faster by analyzing a subset of genes
dge <- dge[1:1000, ]

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~ Disease + (1 | Individual)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, metadata)
fitmm <- eBayes(fitmm)

# Examine design matrix
head(fitmm$design, 3)

# Get results of hypothesis test on coefficients of interest
topTable(fitmm, coef = "Disease1", number = 3)

#Cutoff of adjusted p vaulues, all genes

topTable(fitmm, coef = "Disease1", number = Inf, p.value = 0.05)
savedTable <- topTable(fitmm, coef = "Disease1", number = Inf, p.value = 0.05) %>%
  arrange(P.Value)

#testing an mds plot?
p <- plotMDS(dge)
#it works on the normalized genes
#you need to define a pch and color vector
pch_list <- c(rep(1, 12), rep(2,12))
pch_color <- c(rep("black",12), rep("red",12))
p <- plotMDS(dge, pch = pch_list, col = pch_color)
#this labels them correctly. Now you know how to label by group!
#Make sure that dge is sorted appropriately
