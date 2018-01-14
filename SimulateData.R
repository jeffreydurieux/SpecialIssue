# Thu Jan 11 14:47:02 2018 ------------------------------
# generate data script


###### generate signals S ########

Sim_Signals <- function(n_signals, samples, type = c("b","gauss")){
  
  if(type == "b"){
    res <- replicate(n = n_signals, expr = icasamp("b", query = "rnd", nsamp = samples))  
  }else{
    res <- replicate(n = n_signals, expr = rnorm(samples))
  }
  
  return(res)
}


###### generate time courses A #######

Sim_A <- function(nscan, n_signals){
  TS <- matrix(data = NA, nrow = nscan, ncol = n_signals)
  for(nts in 1:n_signals){# add a piece of code such that NaN's do not occur
    options(warn = -1) #temp warning supression
    repeat{
      TS[ ,nts] <- neuRosim::simTSrestingstate(nscan = nscan, TR = 2, noise = "none")
      if(any(is.finite(TS[,nts])) == TRUE){
        break
      }#end break
    }#end repeat
  }# end for nts
  options(warn = 0) #turn warnings on
  return(TS)
}


########### mix data #############

Sim_X <- function(n_clusters=2, c_size="equal", n_signals = 10, gauss=0, error=0.05, correlation=.7){
  # if gauss
  
  covariance <- matrix(data = correlation, nrow = n_clusters, ncol = n_clusters)
  diag(covariance) <- 1
  
  if(gauss == 1){
    if(n_clusters == 2 & c_size == "equal"){
      A1 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A2 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A <- c(A1,A2)
      
      Al1 <- replicate(30,A[[1]], simplify = FALSE)
      Al2 <- replicate(30,A[[2]], simplify = FALSE)
      Al1 <- lapply(Al1, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al2 <- lapply(Al2, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      
      
      Sa <- clusterwise_correlation(nclusters = 2,covariance = covariance
                                     ,sig = n_signals/2,samples = 2000, type = "b")
      Sb <- clusterwise_correlation(nclusters = 2,covariance = covariance
                                     ,sig = n_signals/2,samples = 2000, type = "gauss")
      
      S1 <- cbind(Sa$SignalList[[1]], Sb$SignalList[[1]])
      S2 <- cbind(Sa$SignalList[[2]], Sb$SignalList[[2]])
      S <- list(S1,S2)
      
      X1 <- lapply(1:length(Al1), function(x) S1 %*% t(Al1[[x]]) )
      X2 <- lapply(1:length(Al2), function(x) S2 %*% t(Al2[[x]]) )
      
      X <- c(X1,X2)
    }else if(n_clusters == 2 & c_size == "unequal"){
      A1 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A2 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A <- c(A1,A2)
      Al1 <- replicate(45,A[[1]], simplify = FALSE)
      Al2 <- replicate(15,A[[2]], simplify = FALSE)
      Al1 <- lapply(Al1, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al2 <- lapply(Al2, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      
      Sa <- clusterwise_correlation(nclusters = 2,covariance = covariance
                                    ,sig = n_signals/2,samples = 2000, type = "b")
      Sb <- clusterwise_correlation(nclusters = 2,covariance = covariance
                                    ,sig = n_signals/2,samples = 2000, type = "gauss")
      
      S1 <- cbind(Sa$SignalList[[1]], Sb$SignalList[[1]])
      S2 <- cbind(Sa$SignalList[[2]], Sb$SignalList[[2]])
      S <- list(S1,S2)
      
      X1 <- lapply(1:length(Al1), function(x) S1 %*% t(Al1[[x]]) )
      X2 <- lapply(1:length(Al2), function(x) S2 %*% t(Al2[[x]]) )
      
      X <- c(X1,X2)
    }else if(n_clusters == 4 & c_size == "equal"){
      A1 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A2 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A3 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A4 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A <- c(A1,A2,A3,A4)
      Al1 <- replicate(15,A[[1]], simplify = FALSE)
      Al2 <- replicate(15,A[[2]], simplify = FALSE)
      Al3 <- replicate(15,A[[3]], simplify = FALSE)
      Al4 <- replicate(15,A[[4]], simplify = FALSE)
      
      Al1 <- lapply(Al1, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al2 <- lapply(Al2, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al3 <- lapply(Al3, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al4 <- lapply(Al4, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      
      
      Sa <- clusterwise_correlation(nclusters = 4,covariance = covariance
                                    ,sig = n_signals/2,samples = 2000, type = "b")
      Sb <- clusterwise_correlation(nclusters = 4,covariance = covariance
                                    ,sig = n_signals/2,samples = 2000, type = "gauss")
      
      S1 <- cbind(Sa$SignalList[[1]], Sb$SignalList[[1]])
      S2 <- cbind(Sa$SignalList[[2]], Sb$SignalList[[2]])
      S3 <- cbind(Sa$SignalList[[3]], Sb$SignalList[[3]])
      S4 <- cbind(Sa$SignalList[[4]], Sb$SignalList[[4]])
      S <- list(S1,S2,S3,S4)
      
      X1 <- lapply(1:length(Al1), function(x) S1 %*% t(Al1[[x]]) )
      X2 <- lapply(1:length(Al2), function(x) S2 %*% t(Al2[[x]]) )
      X3 <- lapply(1:length(Al3), function(x) S3 %*% t(Al3[[x]]) )
      X4 <- lapply(1:length(Al4), function(x) S4 %*% t(Al4[[x]]) )
      
      X <- c(X1,X2,X3,X4)
    }else if(n_clusters ==4 & c_size == "unequal"){
      A1 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A2 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A3 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A4 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A <- c(A1,A2,A3,A4)
      Al1 <- replicate(20,A[[1]], simplify = FALSE)
      Al2 <- replicate(20,A[[2]], simplify = FALSE)
      Al3 <- replicate(10,A[[3]], simplify = FALSE)
      Al4 <- replicate(10,A[[4]], simplify = FALSE)
      Al1 <- lapply(Al1, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al2 <- lapply(Al2, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al3 <- lapply(Al3, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al4 <- lapply(Al4, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      
      
      
      Sa <- clusterwise_correlation(nclusters = 4,covariance = covariance
                                    ,sig = n_signals/2,samples = 2000, type = "b")
      Sb <- clusterwise_correlation(nclusters = 4,covariance = covariance
                                    ,sig = n_signals/2,samples = 2000, type = "gauss")
      
      S1 <- cbind(Sa$SignalList[[1]], Sb$SignalList[[1]])
      S2 <- cbind(Sa$SignalList[[2]], Sb$SignalList[[2]])
      S3 <- cbind(Sa$SignalList[[3]], Sb$SignalList[[3]])
      S4 <- cbind(Sa$SignalList[[4]], Sb$SignalList[[4]])
      S <- list(S1,S2,S3,S4)
      
      
      X1 <- lapply(1:length(Al1), function(x) S1 %*% t(Al1[[x]]) )
      X2 <- lapply(1:length(Al2), function(x) S2 %*% t(Al2[[x]]) )
      X3 <- lapply(1:length(Al3), function(x) S3 %*% t(Al3[[x]]) )
      X4 <- lapply(1:length(Al4), function(x) S4 %*% t(Al4[[x]]) )
      
      X <- c(X1,X2,X3,X4)
    }
    
  }else if(gauss == 0){
    if(n_clusters == 2 & c_size == "equal"){
      A1 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A2 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A <- c(A1,A2)
      Al1 <- replicate(30,A[[1]], simplify = FALSE)
      Al2 <- replicate(30,A[[2]], simplify = FALSE)
      Al1 <- lapply(Al1, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al2 <- lapply(Al2, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      
      
      Sa <- clusterwise_correlation(nclusters = 2,covariance = covariance
                                    ,sig = n_signals,samples = 2000, type = "b")
      
      S1 <- Sa$SignalList[[1]]
      S2 <- Sa$SignalList[[2]]
      S <- list(S1,S2)
      
      X1 <- lapply(1:length(Al1), function(x) S1 %*% t(Al1[[x]]) )
      X2 <- lapply(1:length(Al2), function(x) S2 %*% t(Al2[[x]]) )
      
      X <- c(X1,X2)
      
    }else if(n_clusters == 2 & c_size == "unequal"){
      A1 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A2 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A <- c(A1,A2)
      Al1 <- replicate(45,A[[1]], simplify = FALSE)
      Al2 <- replicate(15,A[[2]], simplify = FALSE)
      Al1 <- lapply(Al1, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al2 <- lapply(Al2, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      
      Sa <- clusterwise_correlation(nclusters = 2,covariance = covariance
                                    ,sig = n_signals,samples = 2000, type = "b")
      
      S1 <- Sa$SignalList[[1]]
      S2 <- Sa$SignalList[[2]]
      S <- list(S1,S2)
      
      X1 <- lapply(1:length(Al1), function(x) S1 %*% t(Al1[[x]]) )
      X2 <- lapply(1:length(Al2), function(x) S2 %*% t(Al2[[x]]) )
      
      X <- c(X1,X2)
      
    }else if(n_clusters == 4 & c_size == "equal"){
      A1 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A2 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A3 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A4 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A <- c(A1,A2,A3,A4)
      Al1 <- replicate(15,A[[1]], simplify = FALSE)
      Al2 <- replicate(15,A[[2]], simplify = FALSE)
      Al3 <- replicate(15,A[[3]], simplify = FALSE)
      Al4 <- replicate(15,A[[4]], simplify = FALSE)
      Al1 <- lapply(Al1, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al2 <- lapply(Al2, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al3 <- lapply(Al3, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al4 <- lapply(Al4, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      
      Sa <- clusterwise_correlation(nclusters = 4,covariance = covariance
                                    ,sig = n_signals,samples = 2000, type = "b")
      
      S1 <- Sa$SignalList[[1]]
      S2 <- Sa$SignalList[[2]]
      S3 <- Sa$SignalList[[3]]
      S4 <- Sa$SignalList[[4]]
      S <- list(S1,S2,S3,S4)
      
      X1 <- lapply(1:length(Al1), function(x) S1 %*% t(Al1[[x]]) )
      X2 <- lapply(1:length(Al2), function(x) S2 %*% t(Al2[[x]]) )
      X3 <- lapply(1:length(Al3), function(x) S3 %*% t(Al3[[x]]) )
      X4 <- lapply(1:length(Al4), function(x) S4 %*% t(Al4[[x]]) )
      
      X <- c(X1,X2,X3,X4)
      
    }else if(n_clusters ==4 & c_size == "unequal"){
      A1 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A2 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A3 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A4 <- replicate(1, expr = Sim_A(nscan = 100, n_signals = n_signals), 
                      simplify = FALSE)
      A <- c(A1,A2,A3,A4)
      Al1 <- replicate(20,A[[1]], simplify = FALSE)
      Al2 <- replicate(20,A[[2]], simplify = FALSE)
      Al3 <- replicate(10,A[[3]], simplify = FALSE)
      Al4 <- replicate(10,A[[4]], simplify = FALSE)
      Al1 <- lapply(Al1, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al2 <- lapply(Al2, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al3 <- lapply(Al3, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      Al4 <- lapply(Al4, FUN = addError, error = 0.5, type = "Gaussian", 
                    additiontype = 1)
      
      Sa <- clusterwise_correlation(nclusters = 4,covariance = covariance
                                    ,sig = n_signals,samples = 2000, type = "b")
      
      S1 <- Sa$SignalList[[1]]
      S2 <- Sa$SignalList[[2]]
      S3 <- Sa$SignalList[[3]]
      S4 <- Sa$SignalList[[4]]
      S <- list(S1,S2,S3,S4)
      
      X1 <- lapply(1:length(Al1), function(x) S1 %*% t(Al1[[x]]) )
      X2 <- lapply(1:length(Al2), function(x) S2 %*% t(Al2[[x]]) )
      X3 <- lapply(1:length(Al3), function(x) S3 %*% t(Al3[[x]]) )
      X4 <- lapply(1:length(Al4), function(x) S4 %*% t(Al4[[x]]) )
      
      X <- c(X1,X2,X3,X4)
    }
    
  }# end if gauss == 0 
  
  # P partitioning vector
  
  if(n_clusters == 2 & c_size == "unequal"){
    P <- c(rep(1,45), rep(2,15))
  }else if(n_clusters == 2 & c_size == "equal"){
    P <- c(rep(1,30), rep(2,30))
  }else if(n_clusters == 4 & c_size == "unequal"){
    P <- c(rep(1,20), rep(2,20), rep(3,10), rep(4,10))
  }else if(n_clusters == 4 & c_size == "equal"){
    P <- c(rep(1,15), rep(2,15), rep(3,15), rep(4,15))
  }
  
  
  Xe <- lapply(X, FUN = addError, error = error, type = "Gaussian", 
              additiontype = 1)
  
  
  res <- list()
  
  if(n_clusters == 4){
    A <-c(Al1,Al2,Al3,Al4)
  }else{
    A <- c(Al1,Al2)
  }
  res$A <- A
  res$S <- S
  res$X <- X
  res$Xe <- Xe
  res$P <- P
  return(res)
}

# test <- Sim_X(n_clusters = 4,c_size = "equal",n_signals = 10,gauss = 0,error = 0.05,correlation = .7)
# 
# diag(congru(test$S[[1]], test$S[[2]]))
# 
# diag(congru(test$S[[1]], test$S[[2]]))
# diag(congru(test$S[[1]], test$S[[3]]))
# diag(congru(test$S[[1]], test$S[[4]]))
# diag(congru(test$S[[2]], test$S[[3]]))
# diag(congru(test$S[[2]], test$S[[4]]))
# diag(congru(test$S[[3]], test$S[[4]]))
# 
