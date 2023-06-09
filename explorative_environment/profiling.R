# Function that computes minus the log-likelihood
negative_log_likelihood_weibull_normal <- function(step1, step2, par,independent = TRUE,
                                                   covarianceMatrix = NULL,  # Specify if the data is assumed independent in the fitting
                                                   eps = 0.001) {
  T <- length(step1)
  # Unpack parameters for the first observation vector
  mu <- par[1:3]
  sigma <- par[4:6]
  # Parameter for the Weibull distribution of the second observation vector
  shape <- par[7:9] 
  scale <- par[10:12]
  # S.T.P 
  Gamma <- diag(3) 
  Gamma[!Gamma] <- par[13:18]  # Fill non-diagonal entries
  Gamma <- Gamma / rowSums(Gamma)  # Divide by row sums
  # Initial distribution
  delta <- c(par[19], par[20], 1)
  delta <- delta / sum(delta)
  all_probs <- matrix(1, nrow = T, ncol = 3)  # Probabilities of observations from the first distribution
  for (i in 1:3) {
    if(independent){
      step_prob1 <- dweibull(step1, shape[i], scale[i])  # Probability density function for the second distribution
      step_prob2 <- dnorm(step2, mean = mu[i], sd = sigma[i])  # Probability density function for the second distribution
      all_probs[, i] <- step_prob1*step_prob2
    }
    else{
      all_probs[, i] <- gDensityRcpp(step1, step2, shape[i], scale[i], mu[i], sigma[i], covarianceMatrix, eps)
    }
  }
  v <- delta * all_probs[1, ]
  llk <- 0
  
  for (t in 2:T) {
    v <- v %*% Gamma * all_probs[t, ]
    llk <- llk + log(sum(v))
    v <- v / sum(v)  # Scaling to avoid numerical problems
  }
  return(-llk)  # Return the negative log-likelihood
}


fit_weibull_normal_hmmRcppDens <- function(step1, # Values in weibull observations
                                           step2, # Values in normal observations
                                           par, # Concatenation of all parameters
                                           Gamma, #Transition probability matrix
                                           delta, # Initial distribution
                                           independent = TRUE,
                                           covarianceMatrix = NULL,  # Specify if the data is assumed independent in the fitting
                                           eps = 0.001)
{ 
  
  # Transform Gamma0 and delta0 to working scale
  wGamma0 <- Gamma0[!diag(3)]
  wDelta0 <- delta0[-3] / delta0[3]
  par0 <- c(par,  wGamma0, wDelta0)
  
  if(!independent){
    ## Use nlminb() to minimise the negative log-likelihood
    mod <- nlminb(start = par0,
                  objective = negative_log_likelihood_weibull_normal,
                  step1 = step1,  # First observation vector
                  step2 = step2,  # Second observation vector
                  lower = c(rep(-Inf, 3), rep(0, 3), rep(0, 6), rep(0, 8)),
                  upper = c(rep(Inf, 3), rep(Inf, 3), rep(Inf, 6), rep(1, 8)),
                  independent = independent,
                  covarianceMatrix = covarianceMatrix,
                  eps = eps)
  }
  else{
    ## Use nlminb() to minimise the negative log-likelihood
    mod <- nlminb(start = par0,
                  objective = negative_log_likelihood_weibull_normal,
                  step1 = step1,  # First observation vector
                  step2 = step2,  # Second observation vector
                  lower = c(rep(-Inf, 3), rep(0, 3), rep(0, 6), rep(0, 8)),
                  upper = c(rep(Inf, 3), rep(Inf, 3), rep(Inf, 6), rep(1, 8)))
  }
  # Normal distribution
  mu_out <- mod$par[1:3]
  sd_out <- mod$par[4:6] 
  # Weibull distribution
  shape_out = mod$par[7:9]
  scale_out = mod$par[10:12]
  # Transition probability matrix
  Gamma_MLE <- diag(3)
  Gamma_MLE[!Gamma_MLE] <- mod$par[13:18] 
  Gamma_MLE <- Gamma_MLE / apply(Gamma_MLE, 1, sum)
  # Initial distribution
  delta_MLE <- c(mod$par[19], mod$par[20], 1)
  delta_MLE <- delta_MLE / sum(delta_MLE)
  return(
    list(mean_normal = mu_out,
         sd_normal = sd_out,
         shape_weibull = shape_out,
         scale_weibull = scale_out,
         Gamma = Gamma_MLE,
         delta = delta_MLE)
  )
}

fit_weibull_normal_hmmRcppDens(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                               delta0, independent = FALSE, covarianceMatrix = covMat)


fit_weibull_normal_hmmFullRcpp <- function(step1, # Values in weibull observations
                                           step2, # Values in normal observations
                                           par, # Concatenation of all parameters
                                           Gamma, #Transition probability matrix
                                           delta, # Initial distribution
                                           covarianceMatrix = NULL,  # Specify if the data is assumed independent in the fitting
                                           eps = 0.001)
{ 
  objective_function <- function(step1, step2, par){
  mu <- par[1:3]
  sigma <- par[4:6]
  # Parameter for the Weibull distribution of the second observation vector
  shape <- par[7:9] 
  scale <- par[10:12]
  # S.T.P 
  Gamma <- diag(3) 
  Gamma[!Gamma] <- par[13:18]  # Fill non-diagonal entries
  Gamma <- Gamma / rowSums(Gamma)  # Divide by row sums
  # Initial distribution
  delta <- c(par[19], par[20], 1)
  delta <- delta / sum(delta)
  
  lik <- negative_log_likelihood_bivariate_weibull_normal_Rcpp(step1, step2, mu, sigma,
                                                               shape, scale, Gamma, covarianceMatrix, delta, eps)
  return(lik)
}
  
  # Transform Gamma0 and delta0 to working scale
  wGamma0 <- Gamma[!diag(3)]
  wDelta0 <- delta[-3] / delta[3]
  par0 <- c(par,  wGamma0, wDelta0)
  
    ## Use nlminb() to minimise the negative log-likelihood
  mod <- nlminb(start = par0,
                objective = objective_function,
                step1 = step1,  # First observation vector
                step2 = step2,  # Second observation vector
                lower = c(rep(-Inf, 3), rep(0, 3), rep(0, 6), rep(0, 8)),
                upper = c(rep(Inf, 3), rep(Inf, 3), rep(Inf, 6), rep(1, 8)))
  # Normal distribution
  mu_out <- mod$par[1:3]
  sd_out <- mod$par[4:6] 
  # Weibull distribution
  shape_out = mod$par[7:9]
  scale_out = mod$par[10:12]
  # Transition probability matrix
  Gamma_MLE <- diag(3)
  Gamma_MLE[!Gamma_MLE] <- mod$par[13:18] 
  Gamma_MLE <- Gamma_MLE / apply(Gamma_MLE, 1, sum)
  # Initial distribution
  delta_MLE <- c(mod$par[19], mod$par[20], 1)
  delta_MLE <- delta_MLE / sum(delta_MLE)
  return(
    list(mean_normal = mu_out,
         sd_normal = sd_out,
         shape_weibull = shape_out,
         scale_weibull = scale_out,
         Gamma = Gamma_MLE,
         delta = delta_MLE)
  )
}

fit_weibull_normal_hmmFullRcpp(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                               delta0,covarianceMatrix = covMat)

a <- microbenchmark(fit_weibull_normal_hmmRcppDens(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                                                   delta0, independent = FALSE, covarianceMatrix = covMat),
                    fit_weibull_normal_hmmFullRcpp(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                                                   delta0,covarianceMatrix = covMat), times = 20)
