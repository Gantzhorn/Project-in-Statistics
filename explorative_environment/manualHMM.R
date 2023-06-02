library(tidyverse)
library(momentuHMM)
library(Rcpp)
library(microbenchmark)
library(Matrix)

Rcpp::sourceCpp("manualHMM.cpp")
# Helper functions for bivariate weibull-normal density
determinantCalculator <- function(Matrix){
  Matrix[1,1]*Matrix[2,2]-Matrix[2,1]*Matrix[1,2]
}

inverseCalculator <- function(Matrix){
  outputMatrix <- Matrix
  outputMatrix[1,1] <- Matrix[2,2]
  outputMatrix[2,2] <- Matrix[1,1]
  outputMatrix[1,2] <- -Matrix[1,2]
  outputMatrix[2,1] <- -Matrix[2,1]
  return(outputMatrix/determinantCalculator(Matrix))
}

g1Inverse <- function(y, mu, sigma){
  (y-mu)/sigma
}

g2Inverse <- function(y, k, lambda, eps){
  qnorm(max(1-exp(-(y/lambda)^k)-eps, eps))
}

gjacobian <- function(y, k, lambda, sigma, eps){
  k*(y/lambda)^k*exp(-(y/lambda)^k)/
    (dnorm(qnorm(1-exp(-(y/lambda)^k)))*y*sigma+eps)
}

quadratic_form <- function(A, x) {
  (t(x) %*% A %*% x)[1]
}

gDensity <- function(x, y, k, lambda, mu, sigma, covarianceMatrix, eps){
  covMatInverse <- inverseCalculator(covarianceMatrix)
  inputvec <- c(g1Inverse(y, mu, sigma), g2Inverse(x, k, lambda, eps))
  result <- 1/(2*pi*sqrt(determinantCalculator(covarianceMatrix)))*
    gjacobian(x, k, lambda, sigma, eps)*
    exp(-1/2*quadratic_form(covMatInverse, inputvec))
  return(result)
}

g2InverseVec <- function(y, k, lambda, eps) {
  updateY <- 1 - exp(-(y / lambda)^k)
  updateY[updateY == 0] <- eps
  updateY[updateY == 1] <- 1-eps
  qnorm(updateY)
}

gjacobianVec <- function(y, k, lambda, sigma, eps) {
  exponent <- (y / lambda)^k
  exponentialPart <- exp(-exponent)
  updateY <- 1 - exponentialPart
  updateY[updateY == 0] <- eps
  updateY[updateY == 1] <- 1-eps
  k * exponent * exponentialPart / (dnorm(qnorm(updateY)) * y * sigma + eps)
}

quadratic_formVec <- function(A, x) {
  apply(x, 1, function(x) t(x) %*% A %*% x)
}

gDensityVec <- function(x, y, k, lambda, mu, sigma, covarianceMatrix, eps) {
  covMatInverse <- inverseCalculator(covarianceMatrix)
  inputMatrix <- cbind(g1Inverse(y, mu, sigma), g2InverseVec(x, k, lambda, eps))
  result <- 1 / (2 * pi * sqrt(determinantCalculator(covarianceMatrix))) *
    gjacobianVec(x, k, lambda, sigma, eps) *
    exp(-1 / 2 * quadratic_formVec(covMatInverse, inputMatrix))
  return(result)
}

## Simulate realisations from a Markov chain

N <- 3 # Number of states
T <- 100 # Number of realisations
delta <- c(1 / 3, 1 / 3, 1 / 3) # Initial distribution
Gamma <- matrix(c(0.9, 0.05, 0.05, 0.05, 0.9, 0.05, 0.05, 0.05, 0.9), ncol = 3) # Transition probability matrix
s <- rep(NA, times = T) # State vector
s[1] <- sample(1:N, size = 1, prob = delta)
for(t in 2:T) {
  s[t] <- sample(1:N, size = 1, prob = Gamma[s[t - 1],])
}

pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## Simulate realisations from an HMM
k <- c(1, 3, 5)
lambda <- c(1, 2, 4)
mu <- c(-1, 0, 1.25) # Means of the state-dependent distributions
sigma <- c(1.5, 0.75, 1) # Standard deviations of the state-dependent distributions
# Simulate independent data
y1 <- rweibull(T, shape = k[s], lambda[s])
y2 <- rnorm(T, mean = mu[s], sd = sigma[s])


covMat <- matrix(c(1,-0.3, -0.3, 1), nrow = 2)
norm_2d <- MASS::mvrnorm(T, mu = c(0,0), Sigma = covMat, empirical = TRUE)

X <- norm_2d[, 1]
Y <- norm_2d[, 2]

X_prime <- qweibull(pnorm(X), shape = k[s], scale = lambda[s])
Y_prime <- Y*sigma[s] + mu[s]
tibble(x = y1, y = y2, s = factor(s)) %>% ggplot(aes(x = x, y = y, col = s)) + geom_point()
tibble(x = X_prime, y = Y_prime, s = factor(s)) %>% ggplot(aes(x = x, y = y, col = s)) + geom_point()

# Initial parameters
mu0 <- mu # Initial values for the mean
sigma0 <- sigma # Initial values for the standard deviation
k0 <- k
lambda0 <- lambda
Gamma0 <- Gamma # Transition probability matrix
delta0 <- delta

fit_weibull_normal_hmm <- function(step1, # Values in weibull observations
                                   step2, # Values in normal observations
                                   par, # Concatenation of all parameters
                                   Gamma, #Transition probability matrix
                                   delta, # Initial distribution
                                   independent = TRUE,
                                   covarianceMatrix = NULL,  # Specify if the data is assumed independent in the fitting
                                   eps = 0.001)
{ 
  
  # Function that computes minus the log-likelihood
  negative_log_likelihood_weibull_normal <- function(step1, step2, par) {
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
        all_probs[, i] <- gDensityVec(step1, step2, shape[i], scale[i], mu[i], sigma[i], covarianceMatrix, eps)
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
  
  # Transform Gamma0 and delta0 to working scale
  wGamma0 <- Gamma0[!diag(3)]
  wDelta0 <- delta0[-3] / delta0[3]
  par0 <- c(par,  wGamma0, wDelta0)
  
  ## Use nlminb() to minimise the negative log-likelihood
  mod <- nlminb(start = par0,
                objective = negative_log_likelihood_weibull_normal,
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

fit_weibull_normal_hmmRcppDens <- function(step1, # Values in weibull observations
                                   step2, # Values in normal observations
                                   par, # Concatenation of all parameters
                                   Gamma, #Transition probability matrix
                                   delta, # Initial distribution
                                   independent = TRUE,
                                   covarianceMatrix = NULL,  # Specify if the data is assumed independent in the fitting
                                   eps = 0.001)
{ 
  
  # Function that computes minus the log-likelihood
  negative_log_likelihood_weibull_normal <- function(step1, step2, par) {
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
  
  # Transform Gamma0 and delta0 to working scale
  wGamma0 <- Gamma0[!diag(3)]
  wDelta0 <- delta0[-3] / delta0[3]
  par0 <- c(par,  wGamma0, wDelta0)
  
  ## Use nlminb() to minimise the negative log-likelihood
  mod <- nlminb(start = par0,
                objective = negative_log_likelihood_weibull_normal,
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
prepSimIndependent <- prepData(tibble(x = y1, y = y2), coordNames = NULL)

# Independent fit
manualMod1 <- fit_weibull_normal_hmm(y1, y2, c(mu0, sigma0, k0, lambda0), Gamma0, delta0) 
momentuMod1 <- fitHMM(data = prepSimIndependent,
                      nbStates = 3,
                      dist = list(x = "weibull",
                                  y = "norm"),
                      Par0 = list(x = c(k0, lambda0),
                                  y = c(mu0, sigma0))
)

# Correlation fit
prepSimCorrelated <- prepData(tibble(x = X_prime, y = Y_prime), coordNames = NULL)
# Independent fit
momentuMod2 <- fitHMM(data = prepSimCorrelated,
                      nbStates = 3,
                      dist = list(x = "weibull",
                                  y = "norm"),
                      Par0 = list(x = c(k0, lambda0),
                                  y = c(mu0, sigma0))
)
manualMod2 <- fit_weibull_normal_hmm(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0, delta0)
# correlated fit
manualMod3 <- fit_weibull_normal_hmm(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                       delta0, independent = FALSE, covarianceMatrix = covMat) 

manualMod3Fast <- fit_weibull_normal_hmmRcppDens(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                                                 delta0, independent = FALSE, covarianceMatrix = covMat)

comp <- microbenchmark(fit_weibull_normal_hmm(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                                     delta0, independent = FALSE, covarianceMatrix = covMat),
               fit_weibull_normal_hmmRcppDens(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                                              delta0, independent = FALSE, covarianceMatrix = covMat), times = 30)


# Plot the densities stratified by the function used
ggplot(df, aes(x = Time, fill = Function)) +
  geom_density(alpha = 0.5) +
  labs(x = "Run Time", y = "Density", fill = "Function") +
  theme_minimal()

hiddenStateSequence_weibullnormal <- function(step1, # Values in weibull observations
                                   step2, # Values in normal observations
                                   par, # Concatenation of all parameters
                                   Gamma, #Transition probability matrix
                                   delta, # Initial distribution
                                   independent = TRUE,
                                   covarianceMatrix = NULL,  # Specify if the data is assumed independent in the fitting
                                   eps = 0.001)
{ 
  
  # Function that computes minus the log-likelihood
  viterbi <- function(step1, step2, par) {
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
    xi <- matrix(0,nrow=T,ncol=3)
    xi[1,] <- v/sum(v)
    for (t in 2:T) {
      v <- apply(xi[t-1,]*Gamma,2,max)*all_probs[t,]
      xi[t,] <- v/sum(v)
    }
    # most probable state sequence
    stSeq <- rep(NA,T)
    stSeq[T] <- which.max(xi[T,])
    for (t in (T-1):1)
      stSeq[t] <- which.max(Gamma[,stSeq[t+1]]*xi[t,])
    return(stSeq)
  }
  
  # Transform Gamma0 and delta0 to working scale
  wGamma0 <- Gamma0[!diag(3)]
  wDelta0 <- delta0[-3] / delta0[3]
  par0 <- c(par,  wGamma0, wDelta0)
  return(viterbi(step1, step2, par0))
}


# independent
hiddenstates1 <- hiddenStateSequence_weibullnormal(y1, y2, c(mu0, sigma0, k0, lambda0), Gamma, delta)
hiddenStatesMomentu1 <- viterbi(momentuMod1)

hiddenStatesMomentu2 <- viterbi(momentuMod2)
hiddenstates2 <- hiddenStateSequence_weibullnormal(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma, delta)
hiddenstates3 <- hiddenStateSequence_weibullnormal(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma, delta, independent = FALSE,
                                       covMat, eps = 0.0001)

