# simulation script for special issue
# Thu Jan 11 14:13:29 2018 ------------------------------
# Extra Wednesday 18-07-2018: do simulation with 80% noise condition. 
# to do this, change for loop iterator: 61:80 from design matrix
# Extra Wednesday 25-07-2018 add high overlap between .07 and .15
# thursday 02-08-2018 last version all factors with additional levels

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)

replication <- args[1]

library("ica")
library("MatrixCorrelation")
library("mclust")
library("apcluster")
library("cluster")
library("neuRosim")

source(file = "/home/durieuxj/CICACODE/R/CICA_datagen.R")
source(file = "/home/durieuxj/SpecialIssue/SimulateData.R")

set.seed(replication)

############ notes ################

#between design

n_clusters <- c(2,4)
c_size <- c("equal", "unequal")
con <- c(.08, .15, .23, .29, .395)
error <- c(0.1, 0.3, 0.6, 0.8)

design <- expand.grid(n_clusters = n_clusters, c_size = c_size,
                      con = con, error = error)

Des <- list()
for (sim in 1:80){
  # generate data
  data <- do.call(SimData, args = as.list(design[sim, ]) )
  
  # single subject ica
  time_icafast <- proc.time()
  icaL <- lapply(data$Xe, FUN = icafast, nc = 20)
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
  SIMMAT <- (RVsM * 0) + (RVsS * 1)
 
  DISSIMMAT <- as.dist(1 - SIMMAT)
  
  # analyze data affinity propagation
 # time_apq0.5 <- proc.time()
 # apq0.5 <- tryCatch( apcluster(s = SIMMAT), warning = function(w)
 #                    return(list(apcluster(s=SIMMAT),w)))
 # time_apq0.5 <- proc.time() - time_apq0.5
  
  #if( is.list(apq0.5) == TRUE ){
  #  apq0.5lab <- 1:60
  #  for(i in 1:length(apq0.5[[1]]@clusters)){
  #    apq0.5lab <- replace(apq0.5lab, apq0.5[[1]]@clusters[[i]], i)  
  #  }
 # }else{
  #  apq0.5lab <- 1:60
  #  for(i in 1:length(apq0.5@clusters)){
  #    apq0.5lab <- replace(apq0.5lab, apq0.5@clusters[[i]], i)  
  #  }
 # }
  
    
  
  
 # time_apq0.1 <- proc.time()
 # apq0.1 <- tryCatch( apcluster(s = SIMMAT, q = 0.1), warning = function(w)
 #   return(list(apcluster(s=SIMMAT, q = 0.1),w)))
 # time_apq0.1 <- proc.time() - time_apq0.1
  
  
 # if( is.list(apq0.1) == TRUE ){
 #   apq0.1lab <- 1:60
  #  for(i in 1:length(apq0.1[[1]]@clusters)){
   #   apq0.1lab <- replace(apq0.1lab, apq0.1[[1]]@clusters[[i]], i)  
   # }
 # }else{
  #  apq0.1lab <- 1:60
   # for(i in 1:length(apq0.1@clusters)){
    #  apq0.1lab <- replace(apq0.1lab, apq0.1@clusters[[i]], i)  
   # }
 # }
  
 time_apqk <- proc.time()
  apqk <- tryCatch( apclusterK(s = SIMMAT, K= design$n_clusters[sim], verbose =F), warning = function(w)
    return(list(apclusterK(s=SIMMAT, K= design$n_clusters[sim],verbose =F),w)))
  time_apqk <- proc.time() - time_apqk
    
  if( is.list(apqk) == TRUE ){
    apqklab <- 1:60
    for(i in 1:length(apqk[[1]]@clusters)){
      apqklab <- replace(apqklab, apqk[[1]]@clusters[[i]], i)  
    }
  }else{
    apqklab <- 1:60
    for(i in 1:length(apqk@clusters)){
      apqklab <- replace(apqklab, apqk@clusters[[i]], i)  
    }
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
 # ARI$ARIap0.5 <- adjustedRandIndex(data$P, apq0.5lab)
 # ARI$ARIap0.1 <- adjustedRandIndex(data$P, apq0.1lab)
  ARI$ARIapk <- adjustedRandIndex(data$P, apqklab)
  ARI$ARIpam <- adjustedRandIndex(data$P, pam$clustering)
  ARI$ARIhclustcomp <- adjustedRandIndex(data$P, hclust_comp_lab)
  ARI$ARIhclustward <- adjustedRandIndex(data$P, hclust_ward_lab)
  
  ARI <- unlist(ARI)
  
  # collect time objects
  # RVM, RVS, ICAS, AP, pam, hclustcomp, hclustward
  TIME <- list()
  #TIME$apq0.1 <- time_apq0.1
  #TIME$apq0.5 <- time_apq0.5
  TIME$apqk <- time_apqk
  TIME$hclustcomp <- time_hclustcomp
  TIME$hclustward <- time_hclustward
  TIME$pam <- time_pam
  TIME$RVM <- time_RVM
  TIME$RVS <- time_RVS
  TIME$ica <- time_icafast
  
  # congru A
  A <- unlist(data$A, recursive = F)
  conA <- lapply(1:length(A), function(x) 
    apply(abs(congru(A[[x]],icaL[[x]]$M)),2,max))
  
  # congru S
  #apply(abs(congru(data$S[[1]],icaL[[1]]$S)),2,max)
  conS <- list()
  for(i in 1:length(data$S)){
    idx <- which(data$P==i)
    conS[[i]] <- lapply(seq_along(idx), function(x) 
      apply(abs(congru(data$S[[i]],icaL[[idx[x]]]$S)),2,max))
  }
  
  conAmean <- sapply(conA, function(x) mean(x))
  conAsd <- sapply(conA, function(x) sd(x))
  
  conSmean <- numeric()
  conSsd <- numeric()
  for(i in 1:length(conS)){
    conSmean[i] <- mean(sapply(conS[[i]], function(x) mean(x)) )
    conSsd[i] <- mean(sapply(conS[[i]], function(x) sd(x)) )
  }
  
  # add to design frame
  Des[[sim]] <- c(design[sim,], ARI, "meanConA" =mean(conAmean), 
                 "meansdConA" = mean(conAsd),
                 "meanConS" = mean(conSmean), "meansdConS" = mean(conAsd)) 
}



# save stuff on export
setwd("/exports/fsw/durieuxj/SpecialIssue/Simulation2/Reduced/")
# Write output, also possible to first save everything in a list object

Results <- do.call(rbind, lapply(Des, data.frame, stringsAsFactors=FALSE))
save(Results, file = paste("Final_Replication_", replication, ".Rdata", sep=""))





