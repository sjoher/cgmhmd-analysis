# Simulation setup
source("Result simulator.R")
library(ggplot2)
library(parallel)
ncores <- detectCores() - 1

nsim <- 25

result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 10, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores = ncores)
}
saveRDS(result, file = "results_1.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 50, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores = ncores)
}
saveRDS(result, file = "results_2.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 100, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores = ncores)
}
saveRDS(result, file = "results_3.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 500, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores = ncores)
}
saveRDS(result, file = "results_4.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 10, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_5.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 50, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_6.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 100, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_7.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 500, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_8.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 10, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_9.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 50, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_10.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 100, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_11.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 500, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_12.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 10, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_13.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 50, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_14.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 100, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_15.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Random", n = 500, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_16.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 10, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_17.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 50, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_18.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 100, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_19.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 500, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_20.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 10, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_21.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 50, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_22.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 100, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_23.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 500, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_24.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 10, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_25.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 50, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_26.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 100, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_27.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 500, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_28.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 10, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_29.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 50, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_30.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 100, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_31.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Cluster", n = 500, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_32.rds") 
rm(result)



result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 10, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_33.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 50, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_34.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 100, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_35.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 500, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_36.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 10, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_37.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 50, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_38.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 100, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_39.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 500, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_40.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 10, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_41.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 50, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_42.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 100, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_43.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 500, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_44.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 10, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_45.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 50, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_46.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 100, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_47.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Approx", network = "Scale-free", n = 500, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_48.rds") 
rm(result)


############################################### Gibbs ###############################################################
#####################################################################################################################
#####################################################################################################################


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 10, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_49.rds") 
rm(result)

#result <- readRDS(file = "results_49.rds")


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 50, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_50.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 100, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_51.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 500, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_52.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 10, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_53.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 50, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_54.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 100, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_55.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 500, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_56.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 10, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_57.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 50, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_58.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 100, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_59.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 500, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_60.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 10, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_61.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 50, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_62.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 100, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_63.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Random", n = 500, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_64.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 10, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_65.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 50, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_66.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 100, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_67.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 500, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_68.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 10, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_69.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 50, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_70.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 100, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_71.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 500, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_72.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 10, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_73.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 50, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_74.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 100, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_75.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 500, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_76.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 10, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_77.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 50, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_78.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 100, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_79.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Cluster", n = 500, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05, nclass = 3, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_80.rds") 
rm(result)



result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 10, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_81.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 50, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_82.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 100, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_83.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 500, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_84.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 10, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_85.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 50, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_86.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 100, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_87.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 500, p = 100, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_88.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 10, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_89.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 50, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_90.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 100, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_91.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 500, p = 50, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_92.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 10, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_93.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 50, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_94.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 100, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_95.rds") 
rm(result)


result <- vector("list", nsim)
for(i in 1:nsim){
  cat(paste(c("This is iteration", i, "out of", length(result), "\n"), collapse =  " "))
  result[[i]] <- simulator(method = "Gibbs", network = "Scale-free", n = 500, p = 100, K = 3, ncat = 6, rho = 1, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, seed = i, ncores=ncores)
}
saveRDS(result, file = "results_96.rds") 
rm(result)
































FPR <- matrix(NA, nrow = nrow(result[[1]]$lambdas), ncol = 6)
TPR <- matrix(NA, nrow = nrow(result[[1]]$lambdas), ncol = 6)
FPR_jgl <- matrix(NA, nrow = nrow(result[[1]]$lambdas), ncol = 6)
TPR_jgl <- matrix(NA, nrow = nrow(result[[1]]$lambdas), ncol = 6)
for(i in 1:6){
  FPR[,i] <- result[[i]]$FPR
  TPR[,i] <- result[[i]]$TPR
  FPR_jgl[,i] <- result[[i]]$FPR_jgl
  TPR_jgl[,i] <- result[[i]]$TPR_jgl
}
fpr_avg <- rowMeans(FPR)
tpr_avg <- rowMeans(TPR)
fpr_avg_jgl <- rowMeans(FPR_jgl)
tpr_avg_jgl <- rowMeans(TPR_jgl)



roc_data <- cbind(fpr_avg, tpr_avg, result[[1]]$lambdas[,2])
roc_data <- as.data.frame(roc_data)
colnames(roc_data) <- c("FPR", "TPR", "l2")
#roc_data$l2 <- as.factor(roc_data$l2)

roc_data_jgl <- cbind(fpr_avg_jgl, tpr_avg_jgl, result[[1]]$lambdas[,2])
roc_data_jgl <- as.data.frame(roc_data_jgl)
colnames(roc_data_jgl) <- c("FPR", "TPR", "l2")
roc_data_jgl[which(roc_data_jgl$l2 == 0),"l2"] <- 0.01
roc_data_jgl[which(roc_data_jgl$l2 == 0.1),"l2"] <- 0.11
roc_data_jgl[which(roc_data_jgl$l2 == 1),"l2"] <- 1.01
roc_data <- rbind(roc_data, roc_data_jgl)
roc_data$l2 <- as.factor(roc_data$l2)
levels(roc_data$l2)
ggplot(data=roc_data,aes(x=FPR, y=TPR, colour=l2)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  scale_color_manual(values=c("#CC0033", "#66CCFF", "#FF0033", "#0099FF", "#FF00CC", "#0000CC"))


fpr_glasso <- matrix(NA, nrow = 6, ncol = length(result[[1]]$FPR_glasso))
tpr_glasso <- matrix(NA, nrow = 6, ncol = length(result[[1]]$FPR_glasso))
for(i in 1:6){
  fpr_glasso[i,] <- result[[i]]$FPR_glasso
  tpr_glasso[i,] <- result[[i]]$TPR_glasso
}
fpr_glasso <- colMeans(fpr_glasso)
tpr_glasso <- colMeans(tpr_glasso)


plot(fpr_glasso, tpr_glasso, type = "l")
