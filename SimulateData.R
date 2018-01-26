# Thu Jan 11 14:47:02 2018 ------------------------------
# generate data script

# con explanation: con .15 is a runif(n, -.15, .15) results in a high correlation between S
# con 1 results in a low correlation

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

SimData <- function(n_clusters, c_size, con, error){
  
  Sbase <- replicate(n = 20, expr = runif(n = 1000, min = -1, max = 1))
  
  Alist <- list()
  Slist <- list()
  Xlist <- list()
  for(i in 1:n_clusters){
    if(c_size == 'equal'){
      Alist[[i]] <- replicate( (60 / n_clusters) , 
                               Sim_A(nscan = 100,n_signals = 20), simplify = F) 
      
      Slist[[i]] <- replicate(n = 20, 
                              runif(n = 1000, min = -con, max = con))
      
    }else if(c_size == 'unequal' & n_clusters == 2){
      Alist[[1]] <- replicate( 15 , Sim_A(nscan = 100,n_signals = 20), simplify = F) 
      Alist[[2]] <- replicate( 45 , Sim_A(nscan = 100,n_signals = 20), simplify = F)
      
      Slist[[i]] <- replicate(n = 20, 
                              runif(n = 1000, min = -con, max = con))
    }else if(c_size == 'unequal' & n_clusters == 4){
      Alist[[1]] <- replicate( 5 , Sim_A(nscan = 100,n_signals = 20), simplify = F) 
      Alist[[2]] <- replicate( 10 , Sim_A(nscan = 100,n_signals = 20), simplify = F) 
      Alist[[3]] <- replicate( 20 , Sim_A(nscan = 100,n_signals = 20), simplify = F) 
      Alist[[4]] <- replicate( 25 , Sim_A(nscan = 100,n_signals = 20), simplify = F) 
      
      Slist[[i]] <- replicate(n = 20, 
                              runif(n = 1000, min = -con, max = con))
    }
    
    S <- lapply(seq_along(Slist), function(x) Sbase + Slist[[x]])
    
    Xlist[[i]] <- lapply(seq_along(Alist[[i]]) , function(x) S[[i]] %*% t(Alist[[i]][[x]]))
    
  }
  
  X <- unlist(Xlist, recursive = F)
  XE <- lapply(X, FUN = addError, error = error)
  
  data <- list()
  data$Xe <- XE
  data$X <- X
  data$S <- S
  data$A <- Alist
  pid <- sapply(Alist, length)
  if(n_clusters == 2){
    data$P <- c(rep(1,pid[1]), rep(2,pid[2]))  
  }else{
    data$P <- c(rep(1,pid[1]), rep(2,pid[2]), rep(3,pid[3]), rep(4,pid[4]))
  }
  
  return(data)
}
