if (!require(here)) install.packages('here')
library(here)

# Load in data
maize <- readRDS("maize_clean.rds")

# Generate sequence of penalty parameters
l1 <- seq(0, 1, 0.1)
l2 <- seq(0, 1, 0.1)

# Estimate networks
est <- cgmmd(maize, "Gibbs", l1, l2)

# Penalty selection 
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
    #lik_mat[i,j] <- (n[j]/2)*(log(det(est[[i]]$Theta[[j]])) - sum(est[[i]]$Theta[[j]] * est[[1]]$Sigma[[j]]))
    lik_mat[i,j] <- n[j]/2 *(determinant(as.matrix(est[[i]]$Theta[[j]]), logarithm = TRUE)$modulus - sum(diag(as.matrix(est[[1]]$ES[[j]]) %*%as.matrix(est[[i]]$Theta[[j]])))) # unpenalized S matrix
    #lik_mat[i,j] <- sum(diag(est[[1]]$Sigma[[j]]*est[[i]]$Theta[[j]])) - log(det(est[[i]]$Theta[[j]]))
    #lik_mat[i,j] <- sum(est[[i]]$Theta[[j]] * est[[1]]$Sigma[[j]]) - log(det(est[[i]]$Theta[[j]]))
    #lik_mat2[i,j] <- n[j]*sum(diag(est[[i]]$Theta[[j]] * est[[1]]$Sigma[[j]])) - n[j]*log(det(est[[i]]$Theta[[j]]))
  }
  lik[i] <- sum(lik_mat[i,])
  aic[i] <- ( - 2 * lik[i] ) + ( 2 * (nu[i,1]+nu[i,2]+nu[i,3]+nu[i,4]) )
  #aic[i] <- ( - 2 * est[[i]]$loglik ) + ( 2 * nu[i,1]+nu[i,2] )
  #bic[i] <- lik_mat[i,1] + nu[i,1]*log(n[1]) + lik_mat[i,2] + nu[i,2]*log(n[2])
  #ebic[i] <- - 2 * est[[i]]$loglik + (log(n[1]) * nu[i,1]) + (log(n[2]) *nu[i,2]) + ( 4 * gamma * log(ncol(maize[[1]])) * nu[i,1]) + (4 * gamma * log(ncol(maize[[1]])) * nu[i,2])
  ebic[i] <- - 2 * lik[i] + (log(n[1]) * nu[i,1]) + (log(n[2]) *nu[i,2]) + 
    (log(n[3]) * nu[i,3]) + (log(n[4]) *nu[i,4]) + 
    (4 * gamma * log(ncol(maize[[1]])) * nu[i,1]) + 
    (4 * gamma * log(ncol(maize[[1]])) * nu[i,2]) +
    (4 * gamma * log(ncol(maize[[1]])) * nu[i,3]) + 
    (4 * gamma * log(ncol(maize[[1]])) * nu[i,4])
}

# AIC penalty
aic_idx <- which.min(abs(aic))
aic_est <- est[[aic_idx]]

# EBIC penalty 
ebic_idx <- which.min(abs(ebic))
ebic_est <- est[[ebic_idx]]

A_aic <- aic_est$Theta
A_ebic <- ebic_est$Theta
G_aic <- vector("list", K)
G_ebic <- vector("list", K)


colnames(maize[[1]])[c(2:6,10,11,13:25,28,29,32:34,37:40,45,54)] <- c("temp_mean", "rad_sum", "evp_sum", "rain_sum", "temp_cv",
                                                                      "max_dry_days", "max_summer_days", 
                                                                         "annual_mean_temp", "rain_wet_q", "rain_dry_q", "rain_warm_q",
                                                                         "rain_cold_q", "temp_range_d", "isothermality",
                                                                         "temp_seas", "max_temp_warm_m", "min_temp_cold_m",
                                                                         "temp_range_y", "mean_temp_wet_q", "mean_temp_dry_q",
                                                                         "cation_ex_cap", "fraction_coarse_frag",
                                                                         "soil_water_ph", "organic_carbon", "organic_carbon_dens",
                                                                         "yield", "yield_theoretical", "growing_degress_",
                                                                         "aridity_idx", "ox_plow_freq", "livestock")
library(igraph)

# Generate graphs using the obtained penalty parameters
for(i in 1:K){
  A_aic[[i]][which(A_aic[[i]] != 0)] <- 1
  A_ebic[[i]][which(A_ebic[[i]] != 0)] <- 1
  colnames(A_aic[[i]]) <- colnames(maize[[1]])
  colnames(A_ebic[[i]]) <- colnames(maize[[1]])
  rownames(A_aic[[i]]) <- colnames(A_aic[[i]])
  rownames(A_ebic[[i]]) <- colnames(A_ebic[[i]])
  colnames(ebic_est$Theta[[i]]) <- colnames(A_ebic[[i]])
  rownames(ebic_est$Theta[[i]]) <- colnames(A_ebic[[i]])
  ebic_est$Theta[[i]] <- ebic_est$Theta[[i]][-c(36,42,55,57),-c(36,42,55,57)]
  A_ebic[[i]] <- A_ebic[[i]][-c(36,42,55,57),-c(36,42,55,57)]
  # remove the unconnected elements from the graphs
  # harvest_doy, stress_pest, soil_conservation, stress_disease
  G_aic[[i]] <- graph_from_adjacency_matrix(A_aic[[i]], mode = "undirected", diag = FALSE)
  G_ebic[[i]] <- graph_from_adjacency_matrix(A_ebic[[i]], mode = "undirected", diag = FALSE)
}

pcor_full <- vector("list", K)
for(k in 1:K){
  pcor_full[[k]] <- ebic_est$Theta[[k]]
}

# Obtain partial correlations
scaled_pcor_full <- pcor_full
for(k in 1:K){
  for(i in 1:nrow(pcor_full[[k]])){
    for(j in 1:ncol(pcor_full[[k]])){
      scaled_pcor_full[[k]][i,j] <- -1*(pcor_full[[k]][i,j]/(sqrt(diag(pcor_full[[k]])[i])*sqrt(diag(pcor_full[[k]])[j])))
    }
  }
  diag(scaled_pcor_full[[k]]) <- diag(pcor_full[[k]])
}

# edge width weights
E(G_ebic[[1]])$e.weight <- 10*abs(scaled_pcor_full[[1]][lower.tri(scaled_pcor_full[[1]])][which(scaled_pcor_full[[1]][lower.tri(scaled_pcor_full[[1]])] !=0)])
E(G_ebic[[2]])$e.weight <- 10*abs(scaled_pcor_full[[2]][lower.tri(scaled_pcor_full[[2]])][which(scaled_pcor_full[[2]][lower.tri(scaled_pcor_full[[2]])] !=0)])
E(G_ebic[[3]])$e.weight <- 10*abs(scaled_pcor_full[[3]][lower.tri(scaled_pcor_full[[3]])][which(scaled_pcor_full[[3]][lower.tri(scaled_pcor_full[[3]])] !=0)])
E(G_ebic[[4]])$e.weight <- 10*abs(scaled_pcor_full[[4]][lower.tri(scaled_pcor_full[[4]])][which(scaled_pcor_full[[4]][lower.tri(scaled_pcor_full[[4]])] !=0)])

# edge color weights
E(G_ebic[[1]])$ec.weight <- scaled_pcor_full[[1]][lower.tri(scaled_pcor_full[[1]])][which(scaled_pcor_full[[1]][lower.tri(scaled_pcor_full[[1]])] !=0)]
E(G_ebic[[2]])$ec.weight <- scaled_pcor_full[[2]][lower.tri(scaled_pcor_full[[2]])][which(scaled_pcor_full[[2]][lower.tri(scaled_pcor_full[[2]])] !=0)]
E(G_ebic[[3]])$ec.weight <- scaled_pcor_full[[3]][lower.tri(scaled_pcor_full[[3]])][which(scaled_pcor_full[[3]][lower.tri(scaled_pcor_full[[3]])] !=0)]
E(G_ebic[[4]])$ec.weight <- scaled_pcor_full[[4]][lower.tri(scaled_pcor_full[[4]])][which(scaled_pcor_full[[4]][lower.tri(scaled_pcor_full[[4]])] !=0)]

# assign color
for(k in 1:K){
  for(i in 1:length(E(G_ebic[[k]])$ec.weight)){
    if(E(G_ebic[[k]])$ec.weight[i] >= 0) E(G_ebic[[k]])$color[i] <- 'green' else E(G_ebic[[k]])$color[i] <- 'red'
  }
}

# Plot all 4 full graphs  
coords <- layout.fruchterman.reingold(G_ebic[[1]])
coordstree <- layout_as_tree(G_ebic[[1]], root = 36)
coords <- rotate_igraph_layout(layout = coordstree, degrees = 270)
plot(G_ebic[[1]],vertex.color="white", vertex.size=degree(G_ebic[[1]])*.35, vertex.label.cex = 0.6, layout = coords, edge.curved = 0.5, edge.color = E(G_ebic[[1]])$color, edge.width = E(G_ebic[[1]])$e.weight, vertex.label.color= "black", vertex.label.dist=1.5, vertex.label.degree=pi)
plot(G_ebic[[2]],vertex.color="white", vertex.size=degree(G_ebic[[2]])*.35, vertex.label.cex = 0.6, layout = coords, edge.curved = 0.5, edge.color = E(G_ebic[[2]])$color, edge.width = E(G_ebic[[2]])$e.weight, vertex.label.color= "black", vertex.label.dist=1.5, vertex.label.degree=pi)
plot(G_ebic[[3]],vertex.color="white", vertex.size=degree(G_ebic[[3]])*.35, vertex.label.cex = 0.6, layout = coords, edge.curved = 0.5, edge.color = E(G_ebic[[3]])$color, edge.width = E(G_ebic[[3]])$e.weight, vertex.label.color= "black", vertex.label.dist=1.5, vertex.label.degree=pi)
plot(G_ebic[[4]],vertex.color="white", vertex.size=degree(G_ebic[[4]])*.35, vertex.label.cex = 0.6, layout = coords, edge.curved = 0.5, edge.color = E(G_ebic[[4]])$color, edge.width = E(G_ebic[[4]])$e.weight, vertex.label.color= "black", vertex.label.dist=1.5, vertex.label.degree=pi)

# Obtain the variables that are connected to yield
linksYactual <- rbind(as.numeric(A_ebic[[1]][,36]), as.numeric(A_ebic[[2]][,36]), 
               as.numeric(A_ebic[[3]][,36]), as.numeric(A_ebic[[4]][,36]))
colnames(linksYactual) <- colnames(A_ebic[[1]])

# Generate yield graph
G_yield <- list(4)
for(i in 1:4){
  G_yield[[i]] <- graph_from_adjacency_matrix(A_ebic[[i]][unique(which(linksYactual == 1, arr.ind = TRUE)[,2]),unique(which(linksYactual == 1, arr.ind = TRUE)[,2])], mode = "undirected", diag = FALSE)
}

coords <- layout.fruchterman.reingold(G_yield[[1]])

# check sign of pcor
# all linked variables
unique(which(linksYactual == 1, arr.ind = TRUE)[,2])
# the colnames of items that should be in the yield graph
cols = colnames(A_ebic[[1]][unique(which(linksYactual == 1, arr.ind = TRUE)[,2]),unique(which(linksYactual == 1, arr.ind = TRUE)[,2])])
pcor <- vector("list", K)
for(k in 1:K){
  pcor[[k]] <- ebic_est$Theta[[k]][cols,cols]
}

# Yield partial correlations
scaled_pcor <- pcor
for(k in 1:K){
  for(i in 1:nrow(pcor[[k]])){
    for(j in 1:ncol(pcor[[k]])){
      scaled_pcor[[k]][i,j] <- -1*(pcor[[k]][i,j]/(sqrt(diag(pcor[[k]])[i])*sqrt(diag(pcor[[k]])[j])))
    }
  }
  diag(scaled_pcor[[k]]) <- diag(pcor[[k]])
}

# assign weight: # 1 is positive 0 is negative
E(G_yield[[1]])$e.weight <- 10*abs(scaled_pcor[[1]][lower.tri(scaled_pcor[[1]])][which(scaled_pcor[[1]][lower.tri(scaled_pcor[[1]])] !=0)])
E(G_yield[[2]])$e.weight <- 10*abs(scaled_pcor[[2]][lower.tri(scaled_pcor[[2]])][which(scaled_pcor[[2]][lower.tri(scaled_pcor[[2]])] !=0)])
E(G_yield[[3]])$e.weight <- 10*abs(scaled_pcor[[3]][lower.tri(scaled_pcor[[3]])][which(scaled_pcor[[3]][lower.tri(scaled_pcor[[3]])] !=0)])
E(G_yield[[4]])$e.weight <- 10*abs(scaled_pcor[[4]][lower.tri(scaled_pcor[[4]])][which(scaled_pcor[[4]][lower.tri(scaled_pcor[[4]])] !=0)])

# edge color weights
E(G_yield[[1]])$ec.weight <- scaled_pcor[[1]][lower.tri(scaled_pcor[[1]])][which(scaled_pcor[[1]][lower.tri(scaled_pcor[[1]])] !=0)]
E(G_yield[[2]])$ec.weight <- scaled_pcor[[2]][lower.tri(scaled_pcor[[2]])][which(scaled_pcor[[2]][lower.tri(scaled_pcor[[2]])] !=0)]
E(G_yield[[3]])$ec.weight <- scaled_pcor[[3]][lower.tri(scaled_pcor[[3]])][which(scaled_pcor[[3]][lower.tri(scaled_pcor[[3]])] !=0)]
E(G_yield[[4]])$ec.weight <- scaled_pcor[[4]][lower.tri(scaled_pcor[[4]])][which(scaled_pcor[[4]][lower.tri(scaled_pcor[[4]])] !=0)]

# assign color
for(k in 1:K){
  for(i in 1:length(E(G_yield[[k]])$ec.weight)){
    if(E(G_yield[[k]])$ec.weight[i] >= 0) E(G_yield[[k]])$color[i] <- 'green' else E(G_yield[[k]])$color[i] <- 'red'
  }
}

# Plot yield graphs
plot(G_yield[[1]], vertex.color="white", vertex.size=degree(G_yield[[1]])+5, vertex.label.cex = 1, layout = coords, edge.color = E(G_yield[[1]])$color, edge.width = E(G_yield[[1]])$e.weight, vertex.label.color= "black", vertex.label.dist=1.2, vertex.label.degree=pi/2)
plot(G_yield[[2]], vertex.color="white", vertex.size=degree(G_yield[[2]])+5, vertex.label.cex = 1, layout = coords, edge.color = E(G_yield[[2]])$color, edge.width = E(G_yield[[2]])$e.weight, vertex.label.color= "black", vertex.label.dist=1.2, vertex.label.degree=pi/2)
plot(G_yield[[3]], vertex.color="white", vertex.size=degree(G_yield[[3]])+5, vertex.label.cex = 1, layout = coords, edge.color = E(G_yield[[3]])$color, edge.width = E(G_yield[[3]])$e.weight, vertex.label.color= "black", vertex.label.dist=1.2, vertex.label.degree=pi/2)
plot(G_yield[[4]], vertex.color="white", vertex.size=degree(G_yield[[4]])+5, vertex.label.cex = 1, layout = coords, edge.color = E(G_yield[[4]])$color, edge.width = E(G_yield[[4]])$e.weight, vertex.label.color= "black", vertex.label.dist=1.2, vertex.label.degree=pi/2)
