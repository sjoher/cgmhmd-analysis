
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


bs1 <- readRDS("bootstrap1.rds")
bs2 <- readRDS("bootstrap2.rds")
bs3 <- readRDS("bootstrap3.rds")
bs4 <- readRDS("bootstrap4.rds")
bs5 <- readRDS("bootstrap5.rds")
bs6 <- readRDS("bootstrap6.rds")
bs7 <- readRDS("bootstrap7.rds")
bs8 <- readRDS("bootstrap8.rds")
bs9 <- readRDS("bootstrap9.rds")
bs10 <- readRDS("bootstrap10.rds")

bsList <- list(bs1, bs2, bs3, bs4, bs5, bs6, bs7, bs8, bs9, bs10)

l1 <- seq(0, 1, 0.1)
l2 <- seq(0, 1, 0.1)

nbootstrap <- 10
ebic_est <- vector("list", nbootstrap)
aic_est <- vector("list", nbootstrap)
for(bs in 1:nbootstrap){
  for(b in 1:20){
    data <- bsList[[bs]][[b]]
    
    n_lambda <- nrow(expand.grid(l1, l2))
    K <- length(maize)
    A1 <- vector("list", n_lambda)
    A2 <- vector("list", n_lambda)
    A3 <- vector("list", n_lambda)
    A4 <- vector("list", n_lambda)
    nu <- matrix(NA, n_lambda, K)
    lik_mat <- matrix(NA, n_lambda, K)
    lik <- numeric(n_lambda)
    aic <- numeric(n_lambda)
    bic <- numeric(n_lambda)
    ebic <- numeric(n_lambda)
    gamma <- 0.5
    n <- sapply(maize, nrow)
    for(i in 1:n_lambda)
    {
      data[[i]]$Theta[[1]] <- as.matrix(data[[i]]$Theta[[1]])
      data[[i]]$Theta[[2]] <- as.matrix(data[[i]]$Theta[[2]])
      data[[i]]$Theta[[3]] <- as.matrix(data[[i]]$Theta[[3]])
      data[[i]]$Theta[[4]] <- as.matrix(data[[i]]$Theta[[4]])
      A1[[i]] <- data[[i]]$Theta[[1]]
      A2[[i]] <- data[[i]]$Theta[[2]]
      A3[[i]] <- data[[i]]$Theta[[3]]
      A4[[i]] <- data[[i]]$Theta[[4]]
      A1[[i]][which(A1[[i]] != 0)] <- 1
      A2[[i]][which(A2[[i]] != 0)] <- 1
      A3[[i]][which(A3[[i]] != 0)] <- 1
      A4[[i]][which(A4[[i]] != 0)] <- 1
      nu[i,1] <- sum(A1[[i]][upper.tri(A1[[i]])])
      nu[i,2] <- sum(A2[[i]][upper.tri(A2[[i]])])
      nu[i,3] <- sum(A3[[i]][upper.tri(A3[[i]])])
      nu[i,4] <- sum(A4[[i]][upper.tri(A4[[i]])])
      for(j in 1:K){
        lik_mat[i,j] <- n[j]/2 *(determinant(as.matrix(data[[i]]$Theta[[j]]), logarithm = TRUE)$modulus - sum(diag(as.matrix(data[[1]]$ES[[j]]) %*%as.matrix(data[[i]]$Theta[[j]])))) # unpenalized S matrix
      }
      lik[i] <- sum(lik_mat[i,])
      aic[i] <- ( - 2 * lik[i] ) + ( 2 * (nu[i,1]+nu[i,2]+nu[i,3]+nu[i,4]) )
      ebic[i] <- - 2 * lik[i] + (log(n[1]) * nu[i,1]) + (log(n[2]) *nu[i,2]) + 
        (log(n[3]) * nu[i,3]) + (log(n[4]) *nu[i,4]) + 
        (4 * gamma * log(ncol(maize[[1]])) * nu[i,1]) + 
        (4 * gamma * log(ncol(maize[[1]])) * nu[i,2]) +
        (4 * gamma * log(ncol(maize[[1]])) * nu[i,3]) + 
        (4 * gamma * log(ncol(maize[[1]])) * nu[i,4])
    }
    
    aic_idx <- which.min(abs(aic))
    aic_est[[bs]][[b]] <- data[[aic_idx]]$Theta
    
    ebic_idx <- which.min(abs(ebic))
    ebic_est[[bs]][[b]] <- data[[ebic_idx]]$Theta
  }
}

ebiclist <- unlist(ebic_est, recursive = FALSE)
length(ebiclist)

# nonzero matrix elements
nonzero <- vector("list", 200)
for(i in 1:200){
  for(k in 1:K){
    nonzero[[i]][[k]] <- which(ebiclist[[i]][[k]][upper.tri(ebiclist[[i]][[k]])] != 0)
  }
}
# we have 200 matrices
# loop over elements of all 200 nonzero matrices
# if encountered, set that element to ++

counter1 <- numeric(1953)
for(i in 1:200){
  for(j in 1:length(nonzero[[i]][[1]])){
      counter1[nonzero[[i]][[1]][j]] <- counter1[nonzero[[i]][[1]][j]] + 1
  }
}
counter2 <- numeric(1953)
for(i in 1:200){
  for(j in 1:length(nonzero[[i]][[2]])){
    counter2[nonzero[[i]][[2]][j]] <- counter2[nonzero[[i]][[2]][j]] + 1
  }
}
counter3 <- numeric(1953)
for(i in 1:200){
  for(j in 1:length(nonzero[[i]][[3]])){
    counter3[nonzero[[i]][[3]][j]] <- counter3[nonzero[[i]][[3]][j]] + 1
  }
}
counter4 <- numeric(1953)
for(i in 1:200){
  for(j in 1:length(nonzero[[i]][[4]])){
    counter4[nonzero[[i]][[4]][j]] <- counter4[nonzero[[i]][[4]][j]] + 1
  }
}

# per edge percentage of how often it is chosen in 200 bootstraps
countList <- list(counter1, counter2, counter3, counter4)
percList <- vector("list", K)
for(k in 1:K){
  for(i in 1:length(countList[[k]])){
    percList[[k]][i] <-countList[[k]][i]/200
  }
}
countList[[1]]
percList[[1]]
which(percList[[1]] > 0.9)


# compare with original results
est <- readRDS("est.rds")

n_lambda <- nrow(expand.grid(l1, l2))
K <- length(maize)
A1 <- vector("list", n_lambda)
A2 <- vector("list", n_lambda)
A3 <- vector("list", n_lambda)
A4 <- vector("list", n_lambda)
nu <- matrix(NA, n_lambda, K)
lik_mat <- matrix(NA, n_lambda, K)
lik <- numeric(n_lambda)
aic <- numeric(n_lambda)
bic <- numeric(n_lambda)
ebic <- numeric(n_lambda)
gamma <- 0.5
n <- sapply(maize, nrow)
for(i in 1:n_lambda)
{
  est[[i]]$Theta[[1]] <- as.matrix(est[[i]]$Theta[[1]])
  est[[i]]$Theta[[2]] <- as.matrix(est[[i]]$Theta[[2]])
  est[[i]]$Theta[[3]] <- as.matrix(est[[i]]$Theta[[3]])
  est[[i]]$Theta[[4]] <- as.matrix(est[[i]]$Theta[[4]])
  A1[[i]] <- est[[i]]$Theta[[1]]
  A2[[i]] <- est[[i]]$Theta[[2]]
  A3[[i]] <- est[[i]]$Theta[[3]]
  A4[[i]] <- est[[i]]$Theta[[4]]
  A1[[i]][which(A1[[i]] != 0)] <- 1
  A2[[i]][which(A2[[i]] != 0)] <- 1
  A3[[i]][which(A3[[i]] != 0)] <- 1
  A4[[i]][which(A4[[i]] != 0)] <- 1
  nu[i,1] <- sum(A1[[i]][upper.tri(A1[[i]])])
  nu[i,2] <- sum(A2[[i]][upper.tri(A2[[i]])])
  nu[i,3] <- sum(A3[[i]][upper.tri(A3[[i]])])
  nu[i,4] <- sum(A4[[i]][upper.tri(A4[[i]])])
  for(j in 1:K){
    lik_mat[i,j] <- n[j]/2 *(determinant(as.matrix(est[[i]]$Theta[[j]]), logarithm = TRUE)$modulus - sum(diag(as.matrix(est[[1]]$ES[[j]]) %*%as.matrix(est[[i]]$Theta[[j]])))) # unpenalized S matrix
  }
  lik[i] <- sum(lik_mat[i,])
  aic[i] <- ( - 2 * lik[i] ) + ( 2 * (nu[i,1]+nu[i,2]+nu[i,3]+nu[i,4]) )
  ebic[i] <- - 2 * lik[i] + (log(n[1]) * nu[i,1]) + (log(n[2]) *nu[i,2]) + 
    (log(n[3]) * nu[i,3]) + (log(n[4]) *nu[i,4]) + 
    (4 * gamma * log(ncol(maize[[1]])) * nu[i,1]) + 
    (4 * gamma * log(ncol(maize[[1]])) * nu[i,2]) +
    (4 * gamma * log(ncol(maize[[1]])) * nu[i,3]) + 
    (4 * gamma * log(ncol(maize[[1]])) * nu[i,4])
}

aic_idx <- which.min(abs(aic))
aic_est <- est[[aic_idx]]

ebic_idx <- which.min(abs(ebic))
ebic_est <- est[[ebic_idx]]

nonzero_orig <- vector("list", K)
for(k in 1:K){
  nonzero_orig[[k]] <- which(ebic_est$Theta[[k]][upper.tri(ebic_est$Theta[[k]])] != 0)
}

# comparing the true positives, i.e. the essential edges, by selecting a cutoff of 0.9
stabtrue <- vector("list", K)
for(k in 1:K){
  for(i in 1:length(which(percList[[k]] > 0.50))){
    stabtrue[[k]][i] <- which(percList[[k]] > 0.50)[i] %in% nonzero_orig[[k]]
  }
}
stabtrue

sum(unlist(stabtrue))/length(unlist(stabtrue))


colnames(ebic_est$Theta[[1]])[c(2:6,10,11,13:25,28,29,32:34,37:40,45,54)] <- c("temp_mean", "rad_sum", "evp_sum", "rain_sum", "temp_cv",
                                                                         "max_dry_days", "max_summer_days", 
                                                                   "annual_mean_temp", "rain_wet_q", "rain_dry_q", "rain_warm_q",
                                                                   "rain_cold_q", "temp_range_d", "isothermality",
                                                                   "temp_seas", "max_temp_warm_m", "min_temp_cold_m",
                                                                   "temp_range_y", "mean_temp_wet_q", "mean_temp_dry_q",
                                                                   "cation_ex_cap", "fraction_coarse_frag",
                                                                   "soil_water_ph", "organic_carbon", "organic_carbon_dens",
                                                                   "yield", "yield_theoretical", "growing_degress_",
                                                                   "aridity_idx", "ox_plow_freq", "livestock")
rownames(ebic_est$Theta[[1]])[c(2:6,10,11,13:25,28,29,32:34,37:40,45,54)] <- c("temp_mean", "rad_sum", "evp_sum", "rain_sum", "temp_cv",
                                                                         "max_dry_days", "max_summer_days", 
                                                                         "annual_mean_temp", "rain_wet_q", "rain_dry_q", "rain_warm_q",
                                                                         "rain_cold_q", "temp_range_d", "isothermality",
                                                                         "temp_seas", "max_temp_warm_m", "min_temp_cold_m",
                                                                         "temp_range_y", "mean_temp_wet_q", "mean_temp_dry_q",
                                                                         "cation_ex_cap", "fraction_coarse_frag",
                                                                         "soil_water_ph", "organic_carbon", "organic_carbon_dens",
                                                                         "yield", "yield_theoretical", "growing_degress_",
                                                                         "aridity_idx", "ox_plow_freq", "livestock")


percmat <- vector("list", k)
for(k in 1:K){
  percmat[[k]] <- matrix(0, 63, 63)
  percmat[[k]][upper.tri(percmat[[k]])][which(ebic_est$Theta[[k]][upper.tri(ebic_est$Theta[[k]])] != 0)] <- percList[[k]][which(ebic_est$Theta[[k]][upper.tri(ebic_est$Theta[[k]])] != 0)]
  colnames(percmat[[k]]) <- colnames(ebic_est$Theta[[1]])
  rownames(percmat[[k]]) <- rownames(ebic_est$Theta[[1]])
  percmat[[k]][lower.tri(percmat[[k]])] = t(percmat[[k]])[lower.tri(percmat[[k]])]
  percmat[[k]] <- percmat[[k]][-c(36,42,55,57),-c(36,42,55,57)]
}

library(reshape2)
melted_percmat <- vector("list", K)
for(k in 1:K){
  melted_percmat[[k]] <- melt(percmat[[k]])
}

library(ggplot2)
ggplot(data = melted_percmat[[1]], aes(x=Var1, y=Var2, fill = value)) +
  scale_fill_gradient2(low = "white", high = "blue", mid = "red", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Bootstrap edge\noccurrence ratio") +
  labs(x = "") +
  labs(y = "") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  geom_tile()

ggplot(data = melted_percmat[[2]], aes(x=Var1, y=Var2, fill = value)) +
  scale_fill_gradient2(low = "white", high = "black", mid = "grey", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Bootstrap edge\noccurrence ratio") +
  labs(x = "") +
  labs(y = "") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  geom_tile()


ggplot(data = melted_percmat[[2]], aes(x=Var1, y=Var2, fill = value)) +
  scale_fill_gradient2(low = "white", high = "blue", mid = "red", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Bootstrap edge\noccurrence ratio") +
  labs(x = "") +
  labs(y = "") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  geom_tile()


ggplot(data = melted_percmat[[3]], aes(x=Var1, y=Var2, fill = value)) +
  scale_fill_gradient2(low = "white", high = "black", mid = "grey", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Bootstrap edge\noccurrence ratio") +
  labs(x = "") +
  labs(y = "") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  geom_tile()


ggplot(data = melted_percmat[[4]], aes(x=Var1, y=Var2, fill = value)) +
  scale_fill_gradient2(low = "white", high = "black", mid = "grey", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Bootstrap edge\noccurrence ratio") +
  labs(x = "") +
  labs(y = "") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  geom_tile()

