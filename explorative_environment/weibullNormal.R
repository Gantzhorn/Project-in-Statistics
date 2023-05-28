library(tidyverse)
M <- 10000

covMat <- matrix(c(1,-0.9, -0.9, 1), nrow = 2)
norm_2d <- MASS::mvrnorm(M, mu = c(0,0), Sigma = covMat, empirical = TRUE)

X <- norm_2d[, 1]
Y <- norm_2d[, 2]

mu_hat <- 4
sigma_hat <- 3

X_prime <- X*sigma_hat + mu_hat

k_hat <- 2
lambda_hat <- 5

Y_prime <- qweibull(pnorm(Y), shape = k_hat, scale = lambda_hat)

# Test of functionality
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

g1 <- function(y, mu, sigma){
  y*mu+sigma
}

g2 <- function(y, k, lambda){
  lambda*(-log(1-pnorm(y)))^(1/k)
}

g1Inverse <- function(y, mu, sigma){
  (y-mu)/sigma
}

g2Inverse <- function(y, k, lambda){
  qnorm(1-exp(-(y/lambda)^k))
}

gjacobian <- function(y, k, lambda, sigma){
  k*(y/lambda)^k*exp(-(y/lambda)^k)/
    (dnorm(qnorm(1-exp(-(y/lambda)^k)))*y*sigma)
}

quadratic_form <- function(A, x) {
  (t(x) %*% A %*% x)[1]
}

gDensity <- function(y, k, lambda, mu, sigma, covarianceMatrix){
  covMatInverse <- inverseCalculator(covarianceMatrix)
  inputvec <- c(g1Inverse(y[1], mu, sigma), g2Inverse(y[2], k, lambda))
  result <- 1/(2*pi*sqrt(determinantCalculator(covarianceMatrix)))*
    gjacobian(y[2], k, lambda, sigma)*
    exp(-1/2*quadratic_form(covMatInverse, inputvec))
  return(result)
}



