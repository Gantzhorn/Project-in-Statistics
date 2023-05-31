library(tidyverse)
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

## Simulate realisations from a Markov chain

N <- 3 # Number of states
T <- 300 # Number of realisations
delta <- c(1 / 3, 1 / 3, 1 / 3) # Initial distribution
Gamma <- matrix(c(0.9, 0.05, 0.05, 0.05, 0.9, 0.05, 0.05, 0.05, 0.9), ncol = 3) # Transition probability matrix
s <- rep(NA, times = T) # State vector
s[1] <- sample(1:N, size = 1, prob = delta)
for(t in 2:T) {
  s[t] <- sample(1:N, size = 1, prob = Gamma[s[t - 1],])
}

pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## Simulate realisations from an HMM





k <- c(1, 7, 20)
lambda <- c(1, 2, 5)
mu <- c(-1, 0, 1) # Means of the state-dependent distributions
sigma <- c(0.5, 0.25, 0.5) # Standard deviations of the state-dependent distributions
# Simulate independent data
y1 <- rweibull(T, shape = k[s], lambda[s])
y2 <- rnorm(T, mean = mu[s], sd = sigma[s])


covMat <- matrix(c(1,-0.9, -0.9, 1), nrow = 2)
norm_2d <- MASS::mvrnorm(T, mu = c(0,0), Sigma = covMat, empirical = TRUE)

X <- norm_2d[, 1]
Y <- norm_2d[, 2]

X_prime <- qweibull(pnorm(X), shape = k[s], scale = lambda[s])
Y_prime <- Y*sigma[s] + mu[s]


# Initial parameters
mu0 <- c(-1, 0, 1) # Initial values for the mean
sigma0 <- c(0.5, 0.25, 0.5) # Initial values for the standard deviation
k0 <- c(1, 7, 20)
lambda0 <- c(1, 2, 5)
Gamma0 <- matrix(c(0.9, 0.05, 0.05, 0.05, 0.9, 0.05, 0.05, 0.05, 0.9), ncol = 3) # Transition probability matrix
delta0 <- c(1 / 3, 1 / 3, 1 / 3)



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
      for(j in 1:T){
        all_probs[j, i] <- gDensity(step1[j], step2[j], shape[i], scale[i], mu[i], sigma[i], covarianceMatrix, eps)
      }
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
              step1 = y1,  # First observation vector
              step2 = y2,  # Second observation vector
              lower = c(rep(-Inf, 3), rep(0, 3), rep(0, 6), rep(0, 8)),
              upper = c(rep(Inf, 3), rep(Inf, 3), rep(Inf, 6), rep(1, 8)),
              control = list(iter.max = 1, eval.max = 1))

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
# Call the function to test
#fit_weibull_normal_hmm(y1, y2, c(mu0, sigma0, k0, lambda0), Gamma0, delta0) # Independent fit

fit_weibull_normal_hmm(X_prime, Y_prime, c(mu0, sigma0, k0, lambda0), Gamma0,
                       delta0, independent = FALSE, covarianceMatrix = covMat)

gDensity(X_prime[2], Y_prime[2], k[1], lambda[1], mu[1], sigma[1], covMat, eps = 0.001)

for(i in 1:3){
  print(i)
  print(gjacobian(X_prime[i], k[2], lambda[2], sigma[2]))
  print(exp(-1/2*quadratic_form(inverseCalculator(covMat), 
                                c(g1Inverse(Y_prime[i],  mu[2], sigma[2]), g2Inverse(X_prime[i], k[2], lambda[2], eps = 0.001)))))
  print(gDensity(X_prime[i], Y_prime[i], k[2], lambda[2], mu[2], sigma[2], covMat, eps = 0.001))
}
