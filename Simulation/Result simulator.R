# Simulator for results
source("Penalty selection.R")
source("Approx method.R")
source("Gibbs method.R")
source("FGL.R")
source("Data simulator.R")
source("JGL.R")

library(MASS)
library(igraph)

simulator <- function(method, network, n, p, K, ncat, rho, gamma_o, gamma_b, gamma_p, prob = NULL, nclass = NULL, seed, ncores){
  set.seed(seed)
  
  sim_dat <- data_sim(network = network, n = n, p = p, K = K, ncat = ncat, rho = rho, gamma_o = gamma_o, gamma_b = gamma_b, gamma_p = gamma_p, prob = prob, nclass = nclass)
  Y <- sim_dat$z
  Theta_0 <- sim_dat$theta
  K <- length(Y)
  lambda1 <- seq(0, 1, 0.05)
  lambda2 <- c(0, 0.1, 1)
  lambdas <- expand.grid(lambda1, lambda2)
  n_lambda <- nrow(lambdas)
  
  if(method == "Gibbs"){
    est = vector("list", n_lambda)
    
    for(combination in 1:n_lambda) 
    {
      if(combination == 1)
      {
        est[[combination]] = vector("list", n_lambda)
        Theta <- vector(mode = "list", length = K)
        for(k in 1:K){
          Theta[[k]] <- diag(ncol(Y[[k]]))
          Theta[[k]] <- as(Theta[[k]], "sparseMatrix") 
        }
        m <- paste(c("Obtaining parameter estimates using the Gibbs method:", floor(100 * combination/n_lambda), "%", "\n"), collapse="")
        cat(m, "\r")
        est[[combination]] = Gibbs_method(Y, lambdas=lambdas, n_lambda=n_lambda, Theta = Theta, combination = combination, max.elongation = 5, em.tol=0.001, ncores = ncores)
      }else{
        est[[combination]] = vector("list", n_lambda)
        for(k in 1:K){
          Theta[[k]] = est[[(combination - 1)]]$Theta[[k]]
          Theta[[k]] = as(Theta[[k]], "dgTMatrix") 
          Theta[[k]] = as(Theta[[k]], "sparseMatrix")
        }    
        m <- paste(c("Obtaining parameter estimates using the Gibbs method:", floor(100 * combination/n_lambda), "%", "\n"), collapse="")
        cat(m, "\r")
        est[[combination]] = Gibbs_method(Y, lambdas=lambdas, n_lambda=n_lambda, Theta= Theta, combination = combination, max.elongation = 5, em.tol=0.001, ncores = ncores)
      }
    }
    for(combination in 1:n_lambda){
      for(k in 1:K){
        est[[combination]]$Theta[[k]] <- as.matrix(est[[combination]]$Theta[[k]])
        est[[combination]]$Sigma[[k]] <- as.matrix(est[[combination]]$Sigma[[k]])
        est[[combination]]$ES[[k]] <- as.matrix(est[[combination]]$ES[[k]])
      }
    }
  } else if(method == "Approx"){
    ini <- initialize(Y, ncores = ncores)
    Z	= ini$Z
    ES	= ini$ES
    lower_upper = ini$lower_upper
    
    est = vector("list", n_lambda)
    for(combination in 1:n_lambda) 
    {
      m <- paste(c("Obtaining parameter estimates using the approximate method:", floor(100 * combination/n_lambda), "%", "\n"), collapse="")
      cat(m, "\r")
      est[[combination]] = vector("list", n_lambda)
      est[[combination]] = Approx_method(Y, Z, ES=ES, lambdas=lambdas, lower_upper=lower_upper, combination = combination, em_tol=0.001, em_iter=5, ncores = ncores)
    }
  } else if(method == "npn"){
    gamma_k <- vector("list", K)
    for(k in 1:K){
      gamma_k[[k]] <- 2*sin(pi/6*cor(Y[[k]], method="spearman"))
    }

    est = vector("list", n_lambda)
    for(combination in 1:n_lambda) 
    {
      m <- paste(c("Obtaining parameter estimates using the npn method:", floor(100 * combination/n_lambda), "%", "\n"), collapse="")
      cat(m, "\r")
      est[[combination]] = vector("list", n_lambda)
      est[[combination]] = fused_gl(S = gamma_k, Y = Y,lambdas[combination,1], lambdas[combination,2])
      est[[combination]]$Theta <- est[[combination]]$theta
    }
  }
  
  # Use the joint graphical lasso method for comparison
  jgl_est <- vector("list", n_lambda)
  for(combination in 1:n_lambda) 
  {
    m <- paste(c("Obtaining parameter estimates using the method by Danaher et al.:", floor(100 * combination/n_lambda), "%", "\n"), collapse="")
    cat(m, "\r")
    jgl_est[[combination]] = vector("list", n_lambda)
    jgl_est[[combination]] = JGL(Y, penalty = "fused", lambdas[combination,1], lambdas[combination,2], return.whole.theta = TRUE)
  }
  
  # Use the glasso for comparison
  library(huge)
  glasso_res <- vector("list", K)
  for(k in 1:K){
    glasso_res[[k]] <- huge.glasso(Y[[k]], lambda = lambda1) 
  }
  
  Theta_true <- Theta_0
  # Cleaning up true theta by setting values close to 0 to 0
  for(k in 1:K){
    Theta_0[[k]][which(abs(Theta_0[[k]]) < 0.0001)] <- 0
  }
  
  # FPR and TPR according to Hao et al. (2018) counting only lower triangle
  FPR_fgl <- vector("list", n_lambda)
  TPR_fgl <- vector("list", n_lambda)
  FPR_jgl <- vector("list", n_lambda)
  TPR_jgl <- vector("list", n_lambda)
  FPR_glasso <- vector("list", n_lambda)
  TPR_glasso <- vector("list", n_lambda)
  for(combination in 1:n_lambda){
    for(k in 1:K){
      FPR_fgl[[combination]][k] <- length(which((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] == 0) & (est[[combination]]$Theta[[k]][lower.tri(est[[combination]]$Theta[[k]], diag = FALSE)] != 0)))/length(which(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] == 0))
      TPR_fgl[[combination]][k] <- length(which((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] != 0) & (est[[combination]]$Theta[[k]][lower.tri(est[[combination]]$Theta[[k]], diag = FALSE)] != 0)))/length(which(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] != 0))
      FPR_jgl[[combination]][k] <- length(which((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] == 0) & (jgl_est[[combination]]$theta[[k]][lower.tri(jgl_est[[combination]]$theta[[k]], diag = FALSE)] != 0)))/length(which(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] == 0))
      TPR_jgl[[combination]][k] <- length(which((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] != 0) & (jgl_est[[combination]]$theta[[k]][lower.tri(jgl_est[[combination]]$theta[[k]], diag = FALSE)] != 0)))/length(which(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] != 0))
    }
    FPR_fgl[[combination]] <- sum(FPR_fgl[[combination]])/K
    TPR_fgl[[combination]] <- sum(TPR_fgl[[combination]])/K
    FPR_jgl[[combination]] <- sum(FPR_jgl[[combination]])/K
    TPR_jgl[[combination]] <- sum(TPR_jgl[[combination]])/K
  }
  FPR_fgl <- unlist(FPR_fgl, use.names = FALSE)
  TPR_fgl <- unlist(TPR_fgl, use.names = FALSE)
  FPR_jgl <- unlist(FPR_jgl, use.names = FALSE)
  TPR_jgl <- unlist(TPR_jgl, use.names = FALSE)
  
  # Compute FPR and TPR for glasso
  FPR_glasso <- vector("list", length(lambda1))
  TPR_glasso <- vector("list", length(lambda1))
  for(combination in 1:length(lambda1)){
    for(k in 1:K){
      FPR_glasso[[combination]][k] <- length(which((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] == 0) & (glasso_res[[k]]$icov[[combination]][lower.tri(glasso_res[[k]]$icov[[combination]], diag = FALSE)] != 0)))/length(which(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] == 0))
      TPR_glasso[[combination]][k] <- length(which((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] != 0) & (glasso_res[[k]]$icov[[combination]][lower.tri(glasso_res[[k]]$icov[[combination]], diag = FALSE)] != 0)))/length(which(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] != 0))
    }
    FPR_glasso[[combination]] <- sum(FPR_glasso[[combination]])/K
    TPR_glasso[[combination]] <- sum(TPR_glasso[[combination]])/K
  }
  FPR_glasso <- unlist(FPR_glasso, use.names = FALSE)
  TPR_glasso <- unlist(TPR_glasso, use.names = FALSE)
  
  
  # Frobenius and entroopy loss for my method and Danaher method
  FL_fgl <- vector("list", n_lambda)
  FL_jgl <- vector("list", n_lambda)
  EL_fgl <- vector("list", n_lambda)
  EL_jgl <- vector("list", n_lambda)
  for(combination in 1:n_lambda){
    for(k in 1:K){
      FL_fgl[[combination]][k] <- sum((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] - est[[combination]]$Theta[[k]][lower.tri(est[[combination]]$Theta[[k]], diag = FALSE)])^2)/sum(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)])
      FL_jgl[[combination]][k] <- sum((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] - jgl_est[[combination]]$theta[[k]][lower.tri(jgl_est[[combination]]$theta[[k]], diag = FALSE)])^2)/sum(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)])
      EL_fgl[[combination]][k] <- sum(diag(solve(Theta_0[[k]])%*%est[[combination]]$Theta[[k]])) - log(det(solve(Theta_0[[k]])%*%est[[combination]]$Theta[[k]])) - p
      EL_jgl[[combination]][k] <- sum(diag(solve(Theta_0[[k]])%*%jgl_est[[combination]]$theta[[k]])) - log(det(solve(Theta_0[[k]])%*%jgl_est[[combination]]$theta[[k]])) - p
    }
    FL_fgl[[combination]] <- sum(FL_fgl[[combination]])/K
    FL_jgl[[combination]] <- sum(FL_jgl[[combination]])/K
    EL_fgl[[combination]] <- sum(EL_fgl[[combination]])/K
    EL_jgl[[combination]] <- sum(EL_jgl[[combination]])/K
    
  }
  FL_fgl <- unlist(FL_fgl, use.names = FALSE)
  FL_jgl <- unlist(FL_jgl, use.names = FALSE)
  EL_fgl <- unlist(EL_fgl, use.names = FALSE)
  EL_jgl <- unlist(EL_jgl, use.names = FALSE)
  
  
  # Frobenius and entropy loss for glasso
  FL_glasso <- vector("list", length(lambda1))
  EL_glasso <- vector("list", length(lambda1))
  for(combination in 1:length(lambda1)){
    for(k in 1:K){
      FL_glasso[[combination]][k] <- sum((Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)] - glasso_res[[k]]$icov[[combination]][lower.tri(glasso_res[[k]]$icov[[combination]], diag = FALSE)])^2)/sum(Theta_0[[k]][lower.tri(Theta_0[[k]], diag = FALSE)])
      EL_glasso[[combination]][k] <- sum(diag(solve(Theta_0[[k]])%*%glasso_res[[k]]$icov[[combination]])) - log(det(solve(Theta_0[[k]])%*%glasso_res[[k]]$icov[[combination]])) - p
    }
    FL_glasso[[combination]] <- sum(FL_glasso[[combination]])/K
    EL_glasso[[combination]] <- sum(EL_glasso[[combination]])/K
  }
  FL_glasso <- unlist(FL_glasso, use.names = FALSE)
  EL_glasso <- unlist(EL_glasso, use.names = FALSE)
  
  
  # Results
  result <- list()
  result$Theta <- vector("list", n_lambda)
  result$jgl_Theta <- vector("list", n_lambda)
  for(combination in 1:n_lambda){
    result$Theta[[combination]] <- est[[combination]]$Theta
    result$jgl_Theta[[combination]] <- jgl_est[[combination]]$theta
  }
  
  result$glasso_Theta <- vector("list", length(lambda1))
  for(combination in 1:length(lambda1)){
    for(k in 1:K){
      result$glasso_Theta[[combination]][k] <- glasso_res[[k]]$icov[[combination]]  
    }
  }
  result$Y <- Y
  result$Theta_0 <- Theta_true
  result$Theta_0_round <- Theta_0
  result$lambdas <- lambdas
  result$FPR <- FPR_fgl
  result$TPR <- TPR_fgl
  result$FPR_jgl <- FPR_jgl
  result$TPR_jgl <- TPR_jgl
  result$FPR_glasso <- FPR_glasso
  result$TPR_glasso <- TPR_glasso
  result$FL <- FL_fgl
  result$FL_jgl <- FL_jgl
  result$FL_glasso <- FL_glasso
  result$EL <- EL_fgl
  result$EL_jgl <- EL_jgl
  result$EL_glasso <- EL_glasso
  result$gamma_o <- gamma_o
  result$gamma_b <- gamma_b
  result$gamma_p <- gamma_p
  result$gamma_g <- 1- gamma_o - gamma_b - gamma_p
  return(result)
}

