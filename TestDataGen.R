library(CICA)
library(ica)
library(MatrixCorrelation)
library(mclust)
library(apcluster)
library(cluster)


# con explanation: con .15 is a runif(n, -.15, .15) results in a high correlation between S
# con 1 results in a low correlation

n_clusters <- c(2,4)
c_size <- c("equal", "unequal")
con <- c(.15, .4, 1)
error <- c(0.1, 0.50, 0.90)

design <- expand.grid(n_clusters = n_clusters, c_size = c_size,
                      con = con, error = error)

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

test <- do.call(SimData, args = as.list(design[4,]))
test$P

icaL <- lapply(data$Xe, FUN = icafast, nc = 20)

comb <- t(combn(1:60, 2))
M <- lapply(icaL, function(anom) anom$M)
S <- lapply(icaL, function(anom) anom$S)

time_RVM <- proc.time()
RVM <- numeric()
for(i in 1:1770){
  RVM[i] <- RV2(M[[ comb[i,1] ]] , M[[ comb[i,2] ]])
}
time_RVM <- proc.time() - time_RVM

time_RVS <- proc.time()
RVS <- numeric()
for(i in 1:1770){
  RVS[i] <- RV2(S[[ comb[i,1] ]] , S[[ comb[i,2] ]])
  #cat(i,"\n")
}
time_RVS <- proc.time() - time_RVS

resM <- cbind(comb,RVM)
resS <- cbind(comb,RVS)

RVsM <- matrix(data = NA, nrow = 60 , ncol = 60)
RVsS <- matrix(data = NA, nrow = 60 , ncol = 60)

for(i in 1:1770){
  val <- resM[i,]
  RVsM[val[1]  , val[2] ] <- val[3]
  
  val2 <- resS[i,]
  RVsS[val2[1]  , val2[2] ] <- val2[3]
  
}

RVsM[lower.tri(RVsM)] = t(RVsM)[lower.tri(RVsM)]
diag(RVsM) <- 1

RVsS[lower.tri(RVsS)] = t(RVsS)[lower.tri(RVsS)]
diag(RVsS) <- 1

SIMMAT <- (RVsM + RVsS) / 2
SIMMAT

DISSIMMAT <- as.dist(1 - SIMMAT)



time_apq0.5 <- proc.time()
apq0.5 <- apcluster(s = SIMMAT)
time_apq0.5 <- proc.time() - time_apq0.5

apq0.5lab <- 1:60
for(i in 1:length(apq0.5@clusters)){
  apq0.5lab <- replace(apq0.5lab, apq0.5@clusters[[i]], i)  
}


time_apq0.1 <- proc.time()
apq0.1 <- apcluster(s = SIMMAT, q = 0.5)
time_apq0.1 <- proc.time() - time_apq0.1

apq0.1lab <- 1:60
for(i in 1:length(apq0.1@clusters)){
  apq0.1lab <- replace(apq0.1lab, apq0.1@clusters[[i]], i)  
}

# pam
time_pam <- proc.time()
pam <- pam(DISSIMMAT, k = 2)
time_pam <- proc.time() - time_pam

#hclust 
time_hclustcomp <- proc.time()
hclust_comp <- hclust(DISSIMMAT, method = "complete")
hclust_comp_lab <- cutree(hclust_comp, k = 2)
time_hclustcomp <- proc.time() - time_hclustcomp

time_hclustward <- proc.time()
hclust_ward <- hclust(DISSIMMAT, method = "ward.D2")
hclust_ward_lab <- cutree(hclust_ward, k = 2)
time_hclustward <- proc.time() - time_hclustward

# Evaluate simulation
ARI <- list()
ARI$ARIap0.5 <- adjustedRandIndex(data$P, apq0.5lab)
ARI$ARIap0.1 <- adjustedRandIndex(data$P, apq0.1lab)
ARI$ARIpam <- adjustedRandIndex(data$P, pam$clustering)
ARI$ARIhclustcomp <- adjustedRandIndex(data$P, hclust_comp_lab)
ARI$ARIhclustward <- adjustedRandIndex(data$P, hclust_ward_lab)

ARI <- unlist(ARI)
ARI