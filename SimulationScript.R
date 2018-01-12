# simulation script for special issue
# Thu Jan 11 14:13:29 2018 ------------------------------


source(file = "SimulateData.R")
source(file = "clusterwise_correlation_v1.R")



############ notes ################

#between design

clusters <- c(2,4)
clustersize <- c("equal", "unequal")
gauss <- c(1,0)
error <- c(0.05, 0.30, 0.50)

design <- expand.grid(clusters=clusters, clustersize = clustersize,
                      gauss = gauss, error = error)
