# Wed Jan 10 11:13:58 2018 ------------------------------
# Jeffrey Durieux
# script for generating correlated signals (clusterwise) 

library(ica)

clusterwise_correlation <- function(nclusters, covariance, sig, samples=1000, type="b"){
  
  if(dim(covariance)[1] != nclusters){
    stop("stop: number of clusters should equal row or col dim of covariance matrix")
  }
  
  R <- chol(covariance)  
  
  # create a signal list, each list element represents a signal.
  # thus length of signalL is the number of signals
  # the first column vectors belong to eachother
  # the second col vectors belong to eachother etc etc
  signalsL <- list()
  for(i in 1:sig){
    
    if(type == "b"){
      signals <- replicate(nclusters, icasamp("b","rnd", samples))  
    }else{
      signals <- replicate(nclusters, rnorm(samples))
    }
    
    
    signalsL[[i]] <- signals %*% R
  }
  
  # combine in large matrix
  signals <- do.call(cbind, signalsL)
  
  # create the right permutation matrix
  seq <- list()
  for(i in 1:nclusters){
    seq[[i]] <- seq(from=i, to = ncol(signals), nclusters)
  }
  seq <- unlist(seq)
  P <- diag(length(seq))[,seq]
  
  mat <- signals %*% P
  
  idx <- split(1:ncol(mat), ceiling(seq_along(1:ncol(mat)) / sig) )
  SignalList <- lapply(idx, function(x) mat[,x])
  
  res <- list()
  res$SignalList <- SignalList
  res$signal <- mat
  return(res)
}


covariance <- matrix(c(1,.7,.7,.7,
                       .7,1.,.7,.7,
                       .7,.7,1,.7,
                       .7,.7,.7,1), nrow = 4)
covariance

test <- clusterwise_correlation(nclusters = 4,covariance = covariance,sig = 5, samples = 1000)

# check correlation
diag(cor(test$SignalList[[1]] , test$SignalList[[2]]) )
diag(cor(test$SignalList[[1]] , test$SignalList[[3]]) )
diag(cor(test$SignalList[[1]] , test$SignalList[[4]]) )

diag(cor(test$SignalList[[2]] , test$SignalList[[3]]) )


library(MatrixCorrelation)

RV2(test$SignalList[[1]], test$SignalList[[2]])
RV2(test$SignalList[[2]], test$SignalList[[3]])
