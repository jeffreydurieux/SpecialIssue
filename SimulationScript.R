# simulation script for special issue
# Thu Jan 11 14:13:29 2018 ------------------------------

library(ica)
library(MatrixCorrelation)
library(mclust)
library(apcluster)
library(cluster)

library(CICA)
source(file = "SimulateData.R")
source(file = "clusterwise_correlation_v1.R")



############ notes ################

#between design

n_clusters <- c(2,4)
c_size <- c("equal", "unequal")
gauss <- c(1,0)
correlation <- c(0.8, 0.3)
error <- c(0.05, 0.50, 0.80)
weight <- c("equal","S","A")
ncomp <- c(10,5,15)

design <- expand.grid(n_clusters=n_clusters, c_size = c_size,
                      gauss = gauss, correlation = correlation,error = error, weight=weight, ncomp=ncomp)


for(sim in 1:24){
  # generate data
  data <- do.call(Sim_X, args = as.list(design[sim,1:5]) )
  
  # single subject ica
  time_icafast <- proc.time()
  icaL <- lapply(data$Xe, FUN = icafast, nc = design$ncomp[sim])
  time_icafast <- proc.time() - time_icafast
  
  # compute RV
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
  
  
  # weighted sim matrix
  if(design$weight[sim] == 'equal'){
    SIMMAT <- (RVsM + RVsS) / 2
  }else if(design$weight[sim] == 'S'){
    SIMMAT <- (RVsM * .1) + (RVsS * .9)
  }else{
    SIMMAT <- (RVsM * .9) + (RVsS * .1)
  }
  
  DISSIMMAT <- as.dist(1 - SIMMAT)
  
  # analyze data affinity propagation
  time_apq0.5 <- proc.time()
  apq0.5 <- apcluster(s = SIMMAT)
  time_apq0.5 <- proc.time() - time_apq0.5
  
  apq0.5lab <- 1:60
  for(i in 1:length(apq0.5@clusters)){
    apq0.5lab <- replace(apq0.5lab, apq0.5@clusters[[i]], i)  
  }
  
  
  time_apq0.1 <- proc.time()
  apq0.1 <- apcluster(s = SIMMAT, q = 0.1)
  time_apq0.1 <- proc.time() - time_apq0.1
  
  apq0.1lab <- 1:60
  for(i in 1:length(apq0.1@clusters)){
    apq0.1lab <- replace(apq0.1lab, apq0.1@clusters[[i]], i)  
  }
  
  # pam
  time_pam <- proc.time()
  pam <- pam(DISSIMMAT, k = design$n_clusters[sim])
  time_pam <- proc.time() - time_pam
  
  #hclust 
  time_hclustcomp <- proc.time()
  hclust_comp <- hclust(DISSIMMAT, method = "complete")
  hclust_comp_lab <- cutree(hclust_comp, k = design$n_clusters[sim])
  time_hclustcomp <- proc.time() - time_hclustcomp
  
  time_hclustward <- proc.time()
  hclust_ward <- hclust(DISSIMMAT, method = "ward.D2")
  hclust_ward_lab <- cutree(hclust_ward, k = design$n_clusters[sim])
  time_hclustward <- proc.time() - time_hclustward
  
  # Evaluate simulation
  ARI <- list()
  ARI$ARIap0.5 <- adjustedRandIndex(data$P, apq0.5lab)
  ARI$ARIap0.1 <- adjustedRandIndex(data$P, apq0.1lab)
  ARI$ARIpam <- adjustedRandIndex(data$P, pam$clustering)
  ARI$ARIhclustcomp <- adjustedRandIndex(data$P, hclust_comp_lab)
  ARI$ARIhclustward <- adjustedRandIndex(data$P, hclust_ward_lab)
  
  ARI <- unlist(ARI)
  
  # collect time objects
  # RVM, RVS, ICAS, AP, pam, hclustcomp, hclustward
  TIME <- list()
  TIME$apq0.1 <- time_apq0.1
  TIME$apq0.5 <- time_apq0.5
  TIME$hclustcomp <- time_hclustcomp
  TIME$hclustward <- time_hclustward
  TIME$pam <- time_pam
  TIME$RVM <- time_RVM
  TIME$RVS <- time_RVS
  TIME$ica <- time_icafast
  
  # congru A
  conA <- lapply(1:length(data$A), function(x) 
    apply(abs(congru(data$A[[x]],icaL[[x]]$M)),2,max))
  
  # congru S
  #apply(abs(congru(data$S[[1]],icaL[[1]]$S)),2,max)
  conS <- list()
  for(i in 1:length(data$S)){
    idx <- which(data$P==i)
    conS[[i]] <- lapply(seq_along(idx), function(x) 
      apply(abs(congru(data$S[[i]],icaL[[idx[x]]]$S)),2,max))
  }
  
  
  
}

# save stuff on export
setwd("/exports/fsw/durieuxj/")
# Write output, also possible to first save everything in a list object

sapply(conA, function(x) mean(x))
sapply(conS[[2]], function(x) mean(x))
