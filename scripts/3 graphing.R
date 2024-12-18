#barplots variances of principal components
#takes normalized counts as an input. use pca_extract as input
pca_barplot <- function(pca_norm_obj){
  barplot(round(pca_norm_obj$sdev^2/sum(pca_norm_obj$sdev^2)*100,2),las=2,
          names.arg=colnames(pca_norm_obj$x),ylab="% Variance explained",
          xlab="PCA principal components")
}