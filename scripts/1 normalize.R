#Function to normalize count matrix
normalize_reads <- function(data){
  # Need to turn this data into a matrix
  count_matrix <- as.matrix(data[,-(1:2)])
  rownames(count_matrix) <- paste0(data[,1],"=", data[,2])
  
  # filter genes by number of counts
  isexpr <- rowSums(cpm(count_matrix) > 0.1) >= 5
  
  # standard usage of limma voom
  dge <- DGEList(count_matrix[isexpr, ])
  dge <- calcNormFactors(dge)
  
  # TODO REMOVE THIS LATER, subsetting dge to make this fast
  # dge <- dge[1:1000, ]
  
  return(dge)
}


# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

