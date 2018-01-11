#clusterwise correlation v2
# Thu Jan 11 13:44:20 2018 ------------------------------

# Tom's advice about how to generate correlated signal matrices

library(ica)
library(CICA)

S <- replicate(n = 4, expr = icasamp("b","rnd",1000))

clusterwise_correlation_v2 <- function(m1,weight){
  E <- addError(datablock = m1, error = 0.05,type = "Gaussian", additiontype = 1)  
  res <- m1 + (weight * E)
  return(res)
}

test <- clusterwise_correlation_v2(S, 2)

congru(S,test)
