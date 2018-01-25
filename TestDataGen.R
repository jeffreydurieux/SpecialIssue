library(CICA)
library(ica)
library(MatrixCorrelation)
library(mclust)
library(apcluster)
library(cluster)
library(distr)

# D <- DExp(rate = 1)
# D1 <- D+1
# D2 <- D+2
# D5 <- D+5
# 
# 
# Sbase <- replicate(n = 10, expr = r(D)(1000))
# 
# SQ1 <- replicate(n = 10, r(D5)(1000))
# SQ2 <- replicate(n = 10, r(D5)(1000))
# 
# S1 <- Sbase + SQ1
# S2 <- Sbase + SQ2

# RV2(S1,S2)
# apply(abs(congru(S1,S2)), 2,max)

conH <- .15
conM <- .3
conL <- 2

Sbase <- replicate(n = 10, expr = runif(n = 1000, min = -1, max = 1))

SQ1 <- replicate(n = 10, runif(n = 1000, min = -conH, max = conH))
SQ2 <- replicate(n = 10, runif(n = 1000, min = -conH, max = conH))
 
S1 <- Sbase + SQ1
S2 <- Sbase + SQ2
 
RV2(S1,S2)
#apply(abs(congru(S1,S2)), 2,max)

A1 <- replicate(10, Sim_A(nscan = 100,n_signals = 10), simplify = F)
A2 <- replicate(10, Sim_A(nscan = 100,n_signals = 10), simplify = F)

X1 <- lapply(seq_along(A1), function(x) S1%*%t(A1[[x]]))
X2 <- lapply(seq_along(A2), function(x) S2%*%t(A2[[x]]))

X <- c(X1,X2)
Xe <- lapply(X, FUN = addError, error=0.1)

data <- list()
data$Xe <- Xe
data$P <- c(rep(1,10), rep(2,10))

icaL <- lapply(data$Xe, FUN = icafast, nc = 10)

comb <- t(combn(1:20, 2))
M <- lapply(icaL, function(anom) anom$M)
S <- lapply(icaL, function(anom) anom$S)

time_RVM <- proc.time()
RVM <- numeric()
for(i in 1:190){
  RVM[i] <- RV2(M[[ comb[i,1] ]] , M[[ comb[i,2] ]])
}
time_RVM <- proc.time() - time_RVM

time_RVS <- proc.time()
RVS <- numeric()
for(i in 1:190){
  RVS[i] <- RV2(S[[ comb[i,1] ]] , S[[ comb[i,2] ]])
  #cat(i,"\n")
}
time_RVS <- proc.time() - time_RVS

resM <- cbind(comb,RVM)
resS <- cbind(comb,RVS)

RVsM <- matrix(data = NA, nrow = 20 , ncol = 20)
RVsS <- matrix(data = NA, nrow = 20 , ncol = 20)

for(i in 1:190){
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

apq0.5lab <- 1:10
for(i in 1:length(apq0.5@clusters)){
  apq0.5lab <- replace(apq0.5lab, apq0.5@clusters[[i]], i)  
}


time_apq0.1 <- proc.time()
apq0.1 <- apcluster(s = SIMMAT, q = 0.1)
time_apq0.1 <- proc.time() - time_apq0.1

apq0.1lab <- 1:10
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