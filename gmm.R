# use faithful data set
library(mvtnorm) # needed to compute the multivariate densities

X <- scale(faithful) # Example Dataset

gmm <- function(X = NULL, K = 2, seed = NULL, max_iter = 100, change = 1e-6){
  # Initial Values
  X <- as.matrix(X)
  N <- nrow(X)
  mu <- list()
  Sigmas = list()
  set.seed(seed)
  init_s <- round(N/K, 0) # initial sample size for each k
  s_id <- 1:N
  for(k in 1:K){ # initial partitioned data into K clusters
    if(k<K){
      sk_id <- sample(x = s_id, size = init_s, replace = FALSE)
      Xk <- X[sk_id, ]
      mu[[k]] <- colMeans(Xk)
      Sigmas[[k]] <- cov(Xk)
      s_id <- s_id[!(s_id%in%sk_id)]
    }
    else{
      Xk <- X[s_id, ]
      mu[[k]] <- colMeans(Xk)
      Sigmas[[k]] <- cov(Xk)
    }
    #print(nrow(Xk))
  }
  
  pi_K <- 1/K # initial values for prior pi
  lik_old <- -1e6
  lik_list <- numeric()
  Norm_d <- matrix(NA, nrow = N, ncol = K)
  i <- 1
  while(i <= max_iter){
    # E step
    for (k in 1:K) { # compute normal density matrix
      Norm_d[,k] <- dmvnorm(x = X, mean = mu[[k]], sigma = Sigmas[[k]])
    }
    
    pi_times_norm <- sweep(x = Norm_d, MARGIN = 2, # pi*Normal(muk, sigmak)
                           FUN = "*", STATS = pi_K)
    gam_denom <- rowSums(pi_times_norm)
    
    gamma_z <- pi_times_norm/gam_denom # Equation 9.23
    
    N_K <- round(colSums(gamma_z), 0) # Equation 9.27
    
    
    # M step
    
    for (k in 1:K) {
      mu[[k]] <- (1/N_K[k])*colSums(gamma_z[,k]*X) # Equation 9.24
      X_mu <- sweep(x = X, MARGIN = 2, FUN = "-", STATS = mu[[k]])
      Sigmas[[k]] <- (1/N_K[k])*t(gamma_z[,k]*X_mu)%*%X_mu # Equation 9.25
    }
    pi_K <- N_K/N # Equation 9.26
    
    # Compute likelihood
    for (k in 1:K) {
      Norm_d[,k] <- dmvnorm(x = X, mean = mu[[k]], sigma = Sigmas[[k]])
    }
    
    pi_times_norm <- sweep(x = Norm_d, MARGIN = 2, FUN = "*", STATS = pi_K)
    lik <- sum(log(rowSums(pi_times_norm)))
    lik_list[i] <- lik
    changes <- abs(lik - lik_old)/ abs(lik_old)
    if(changes <= change){
      break
    }
    lik_old <- lik
    i <- i + 1
  }
  out <- list(lik_list = lik_list, 
              gamma_z = gamma_z, mu = mu, 
              Sigmas = Sigmas, pi_K = pi_K, n_iter = i, N_K = N_K)
}

gmm_out <- gmm(X = X, K = 2, seed = 1424, max_iter = 1000)

