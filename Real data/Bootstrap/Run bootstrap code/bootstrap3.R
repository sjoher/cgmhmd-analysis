if (!require(here)) install.packages('here')
library(here)

# this is a boostrap analysis of the maize data
# we resample the rows of the data with replacement to generate 200
# "new" datasets whereby we fit a copula graphical model on the original data
# perform model selection using AIC and EBIC and do the same on the bootstrapped data, select
# a graph using the penalty parameters of the bootstrap iteration
# then we evaluate how often (percentage) certain edges are and are not present in the graphs
# of the bootstrapped graphs and evaluate the discrepancy of these

# mention in data application that we use ebic and aic because they are computationally efficient

# we only use the Gibbs method in this case

maize <- readRDS("maize_clean.rds")

# load copula function
source("cgmmd_boostrap.R")

l1 <- seq(0, 1, 0.1)
l2 <- seq(0, 1, 0.1)
expand.grid(l1, l2)

# we will choose for the nonpar version

est <- vector("list", 20)
for(i in 1:20){
  set.seed(i+40)
  # permute rows
  idx1 <- sample(1:nrow(maize[[1]]), nrow(maize[[1]]), replace = TRUE)
  idx2 <- sample(1:nrow(maize[[2]]), nrow(maize[[2]]), replace = TRUE)
  idx3 <- sample(1:nrow(maize[[3]]), nrow(maize[[3]]), replace = TRUE)
  idx4 <- sample(1:nrow(maize[[4]]), nrow(maize[[4]]), replace = TRUE)
  
  maize[[1]] <- maize[[1]][idx1,]
  maize[[2]] <- maize[[2]][idx2,]
  maize[[3]] <- maize[[3]][idx3,]
  maize[[4]] <- maize[[4]][idx4,]
  
  est[[i]] <- cgmmd(maize, "Gibbs", l1, l2, ncores)
}

saveRDS(est, here("Results", "bootstrap3.rds")) 

