# simulation script for special issue
# Thu Jan 11 14:13:29 2018 ------------------------------

library(ica)
library(CICA)
library(mclust)
source(file = "SimulateData.R")
source(file = "clusterwise_correlation_v1.R")



############ notes ################

#between design

n_clusters <- c(2,4)
c_size <- c("equal", "unequal")
gauss <- c(1,0)
error <- c(0.05, 0.30, 0.50)

design <- expand.grid(n_clusters=n_clusters, c_size = c_size,
                      gauss = gauss, error = error)


for(i in 1:24){
  # generate data
  data <- do.call(Sim_X, args = as.list(design[i,]) )
  
  
  # analyze data
  tmp <- proc.time()
  cica <- CICA(nStarts = 30, DataList = data$Xe, nComp = 10, 
               nClus = design$n_clusters[i], show = "best")
  time <- proc.time() - tmp
  
  
  # Evaluate simulation
  ARI <- adjustedRandIndex(data$P, cica$P)
  
  #### congru of S
  
  #### congru of A
  
  
  # save stuff on export
  setwd("/exports/fsw/durieuxj/")
  
  # Write output, also possible to first save everything in a list object
  save(data, file =paste("Simdata","Row", i, "Rep",Replication ,".Rdata" , sep ="")
  save(Analysis, file =paste("Analysis","Row", i, "Rep",Replication ,".Rdata" , sep ="")
  save(ARI, file =paste("Evaluation","Row", i, "Rep",Replication ,".Rdata" , sep ="")
  save(time, file =paste("Time","Row",i, "Rep",Replication ,".Rdata" , sep ="")
  
}


