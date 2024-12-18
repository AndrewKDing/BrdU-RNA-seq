#extracts principal components from list of DGEs. To be used for plotting
pca_extract <- function(x){

  pca_norm <- cpm(x, log=TRUE, prior.count=3)
  pca_norm_obj <- prcomp(t(pca_norm))
  
  return(pca_norm_obj)
}
