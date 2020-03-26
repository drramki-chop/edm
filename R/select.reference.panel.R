#' @export
select.reference.panel <- function(test.counts, reference.counts, bin.length = NULL, data = NULL, correlations=NULL) {
  
  message('Optimization of the choice of aggregate reference set, this process can take some time')
  
  if (sum(test.counts > 2) < 5) {
    message('It looks like the test samples has only ', sum(test.counts > 2), ' bins with 2 or more reads. The coverage is too small to perform any meaningful inference so no likelihood will be computed.')
    #my.res <- list(reference.choice = dimnames(reference.counts)[[2]][1], summary.stats = NULL)
    #return(my.res)
    
    # we'll return NULL as that's what we check in call variants
    return(NULL)
  }
  
  if (class(reference.counts) != 'matrix') stop('The reference sequence count data must be provided as a matrix')
  if (nrow(reference.counts) != length(test.counts)) stop("The number of rows of the reference matrix must match the length of the test count data\n")
  if (is.null(bin.length)) bin.length <- rep(1, length(test.counts))
  
  reference.counts <- reference.counts[, order(correlations, decreasing = TRUE), drop = FALSE]
  top.columns <- dim(reference.counts)[2]
  if (top.columns > 50){
     top.columns = 50
    }
  csum <- t(apply(reference.counts[,1:top.columns],1,cumsum)) +1
  ratio <- apply(csum,2,function(x){return(mean(test.counts/x,na.rm=T))})
  if( length(which(ratio <0.05)) > 0) {
  my.res <- list(reference.choice = names(ratio)[1:min(which(ratio < 0.05))])
  return(my.res)
  } else {
    message('It looks like the test sample is an outlier and CNV cannot be called.')
    #my.res <- list(reference.choice = dimnames(reference.counts)[[2]][1], summary.stats = NULL)
    return(NULL)
  }
  
}
