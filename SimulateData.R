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

