#let's make some PCA graphs
# normalize.R normalizes differentially expressed genes, we need to take the 
# list of dges to make a pca plot.

#extracts principal components from list of DGEs. To be used for plotting
pca_extract <- function(x){

  pca_norm <- cpm(x, log=TRUE, prior.count=3)
  pca_norm_obj <- prcomp(t(pca_norm))
  
  #pulling out the principal components
  pca_matrix <- as.data.frame(pca_norm_obj$x) 
  
  return(pca_matrix)
}
