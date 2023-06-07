library(tidyverse)
library(ggExtra)
library(xtable)
library(momentuHMM)
theme_set(theme_bw())
library(gridExtra)
plotPath = "/home/anders/Desktop/Skole/Blok 4 - 2023/Project in Statistics/Project-in-Statistics/latex/figures"
proj_palette <- c("#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00",
                  "#CC79A7")
M <- 2500
rhos <- c(-0.9, -0.45, -0.1, 0.9, 0.45, 0.1)
mu_hat <- 4
sigma_hat <- 3
k_hat <- 2
lambda_hat <- 5
cor_vec <- numeric(length(rhos))
k_est <- numeric(length(rhos))
lambda_est <- numeric(length(rhos))

mu_est <- numeric(length(rhos))
sigma_est <- numeric(length(rhos))

plots <- list()

for (i in seq_along(rhos)) {
  rho1 <- rhos[i]

  covMat <- matrix(c(1, rho1, rho1, 1), nrow = 2)
  norm_2d <- MASS::mvrnorm(M, mu = c(0, 0), Sigma = covMat, empirical = TRUE)
  
  X <- norm_2d[, 1]
  Y <- norm_2d[, 2]
  
  X_prime <- X * sigma_hat + mu_hat

  Y_prime <- qweibull(pnorm(Y), shape = k_hat, scale = lambda_hat)
  
  k_est[i] <- uniroot(function(k) {sum(Y_prime^k * log(Y_prime)) / sum(Y_prime^k) - 1/k - mean(log(Y_prime))},
                   c(0.1, 10))$root
  lambda_est[i] <- (mean(Y_prime^k_est[i]))^(1/k_est[i])
  
  mu_est[i] <- mean(X_prime)
  
  sigma_est[i] <- sd(X_prime)
  cor_vec[i] <- cor(X_prime, Y_prime)
  p <- ggplot(tibble(x = Y_prime, y = X_prime), aes(x = x, y = y)) +
    geom_point(alpha = 0.1, color = proj_palette[1]) +
    labs(x = if(rho1>0){expression(bold(Z1))} else{""}, y = if(abs(rho1)==0.9){(expression(bold(Z2)))} else{""}) +
    xlim(0, 17) +
    ylim(-10, 15) +
    theme(axis.title = element_text(face = "bold"))
  
  
  plots[[as.character(rho1)]] <- ggMarginal(p, fill = proj_palette[2])
}

# Create a data frame with the values
parameter_results <- data.frame(rhos, k_est, lambda_est, mu_est, sigma_est, cor_vec)

# Convert the data frame to an xtable object
xtable_table <- xtable(parameter_results, digits = 6)

ggsave("correlatedWeibullNormalSim.jpeg", plot = grid.arrange(grobs = plots, ncol = 3), path = plotPath)

# Actual data fitting
GoldenEagle_FinalData_Dryad_FE <- read_csv(
  "~/Desktop/Skole/Blok 4 - 2023/Project in Statistics/Project-in-Statistics/rmarkdown/GoldenEagle_FinalData_Dryad_FE.csv")

eagleData <- GoldenEagle_FinalData_Dryad_FE %>% group_by(Segment_ID) %>%
  mutate(vertical_steps = Altitude - lag(Altitude, default = first(Altitude))) %>% 
  rename(horizontal_steps = Step_length) %>% 
  slice(2:n()) %>% ungroup()

# The idea:
k_hat <- uniroot(function(k) {sum(eagleData$horizontal_steps^k * log(eagleData$horizontal_steps)) / sum(eagleData$horizontal_steps^k) - 1/k - mean(log(eagleData$horizontal_steps))},
                 c(0.1, 10))$root
lambda_hat <- (mean(eagleData$horizontal_steps^k_hat))^(1/k_hat)

mu_hat <- mean(eagleData$vertical_steps)
sigmasq_hat <- var(eagleData$vertical_steps)

# covMatEagle <- matrix(c(1, -0.3,
#                         -0.3, 1), ncol = 2,
#                       byrow = TRUE)
# 
# 
# bivariateNormEagle = MASS::mvrnorm(n=50000, mu=c(0, 0),
#                                    Sigma=covMatEagle, empirical=TRUE)
# X = bivariateNormEagle[, 1]
# Y = bivariateNormEagle[, 2]
# 
# weibullNormalEagle <- tibble(weibull = qweibull(pnorm(bivariateNormEagle[,2], 0, sqrt(covMatEagle[2,2])),
#                                                 shape = k_hat, scale = lambda_hat),
#                              normal = bivariateNormEagle[,1]*sqrt(sigmasq_hat)+mu_hat)
# 
# 
# ggplot() + geom_point(data = eagleData, mapping = aes(x = horizontal_steps, y = vertical_steps), 
#                       alpha = .1) +
#   stat_density_2d(data = weibullNormalEagle, mapping = aes(x=weibull, y=normal, fill = ..level..)
#                   , geom = "polygon", colour="white") +
#   scale_fill_gradient(low = "darkred", high = "firebrick", na.value = NA)
# 
# # look at marginals
# # Check if the follow the distributions they should according to the simulations
# ggplot() +
#   geom_histogram(data = weibullNormalEagle, mapping = aes(x = weibull, after_stat(density)), 
#                  col = "black", fill = "firebrick") +
#   stat_function(fun = dweibull, args = list(shape = k_hat1, scale = lambda_hat1), col = "black")
# 
# ggplot() +
#   geom_histogram(data = weibullNormalEagle, mapping = aes(x = normal, after_stat(density)), 
#                  col = "black", fill = "firebrick") +
#   stat_function(fun = dnorm, args = list(mu_hat1, sigma_hat1))
# 
# # Check if they follow the actual data
# tibble(data = c(eagleData$horizontal_steps, weibullNormalEagle$weibull), 
#        label = factor(rep(c("Observed", "Simulated"),
#                           times = c(length(eagleData$horizontal_steps), length(weibullNormalEagle$weibull))))
#        ) %>% ggplot(aes(x = data, after_stat(density), col = label)) + geom_density() +
#   scale_fill_manual(values = c("firebrick", "dodgerblue3"), labels = c("Observed Data", "Simulated Data"),
#                     guide = "legend", name = "")
# 
# tibble(data = c(eagleData$vertical_steps, weibullNormalEagle$normal), 
#        label = factor(rep(c("Observed", "Simulated"),
#                           times = c(length(eagleData$vertical_steps), length(weibullNormalEagle$normal))))
# ) %>% ggplot(aes(x = data, after_stat(density), col = label)) + geom_density() +
#   scale_fill_manual(values = c("firebrick", "dodgerblue3"), labels = c("Observed Data", "Simulated Data"),
#                     guide = "legend", name = "")
# Do it based on states
prepEagleData <- momentuHMM::prepData(data = data.frame(ID = eagleData$Segment_ID,
                                                        horizontal_steps = eagleData$horizontal_steps,
                                                        vertical_steps = eagleData$vertical_steps),
                                      coordNames = NULL)

N <- 3 
# Initial values for weibull distributions of horizontal distances
k0 <- rep(k_hat, times = N)
lambda0 <- rep(lambda_hat, times = N)

# Initial values for normal distributions for vertical distances
mu0 <- rep(mu_hat, times = N)
sigma0 <- rep(sigmasq_hat, times = N)
eagleFit1 <- momentuHMM::fitHMM(prepEagleData,
                                nbStates = N,
                                dist = list(horizontal_steps = "weibull",
                                            vertical_steps = "norm"),
                                Par0 = list(horizontal_steps = c(k0, lambda0),
                                            vertical_steps = c(mu0, sigma0))
)
decodedStates <- momentuHMM::viterbi(eagleFit1)
localDecoding <- max.col(momentuHMM::stateProbs(eagleFit1))

obsPerAnimal <- as.numeric(table(prepEagleData$ID))

statesForSim <- rep(0, nrow(eagleFit1$data))
k <- 0
for (i in seq_along(obsPerAnimal)){
  statesForSim[k+1] <- sample.int(N, 1, prob = eagleFit1$mle$delta[1,])
  for(j in 2:obsPerAnimal[i]){
    statesForSim[j+k] <- sample.int(N, 1, prob = eagleFit1$mle$gamma[statesForSim[j+k-1], ])
  }
  k <- obsPerAnimal[i]+k
}
statesForSim[c(1, 1+cumsum(obsPerAnimal)[-length(obsPerAnimal)])]

statesForSim %>% table()
decodedStates %>% table()

weibullParameters <- matrix(c(eagleFit1$mle[[1]][1, ],
                              eagleFit1$mle[[1]][2, ]), ncol = N, byrow = TRUE)
normalParameters <- eagleFit1$mle[[2]]


eagleGamma <- eagleFit1$mle$gamma

# Fit similar data
covMatEagle <- matrix(c(1, 0.99,
                        0.99, 1), ncol = 2,
                      byrow = TRUE)


bivariateNormEagle = MASS::mvrnorm(n=nrow(eagleFit1$data), mu=c(0, 0),
                                   Sigma=covMatEagle, empirical=TRUE)
X = bivariateNormEagle[, 1]
Y = bivariateNormEagle[, 2]

weibullNormalEagle <- tibble(weibull = qweibull(pnorm(X, 0, sqrt(covMatEagle[2,2])),
                                                shape = weibullParameters[1, ][statesForSim],
                                                scale = weibullParameters[2, ][statesForSim]),
                             normal = Y*normalParameters[2, ][statesForSim]+normalParameters[1, ][statesForSim],
                             state = factor(statesForSim))

probStates <- table(statesForSim)/nrow(eagleFit1$data)

x <- density(eagleData$horizontal_steps)$x

pdf_values <- sapply(1:N, function(i) {
  # Take care of scenarios where there are 1 and 2 parameters in distribution
  density <- dweibull(x, weibullParameters[1, ][i], weibullParameters[2, ][i])*probStates[i]
  density <- density
  density
})

df <- data.frame(x = rep(x, (N+1)),
                 density = as.vector(cbind(pdf_values, rowSums(pdf_values))),
                 state = factor(rep(c(1:N, "Total"), each = length(x))))

# Plot the density functions using ggplot2
weightedPlot <- ggplot() + 
  geom_histogram(data = tibble(x = eagleData$horizontal_steps),
                 aes(x = x, after_stat(density)), fill = "gray", col = "black") +
  geom_line(data = filter(df, state != "Total"), aes(x = x, y = density, color = state),
            lwd = 1, alpha = 1.5) +
  geom_line(data = filter(df, state == "Total"), aes(x = x, y = density),
            lwd = 1, alpha = 1.5, linetype = "dashed") +
  xlab("horizontal_steps") +
  ylab("Density")

# Marginally the data looks like this:
ggplot() + geom_density(eagleData, mapping = aes(x = horizontal_steps),
                        col = "dodgerblue4", lwd = 2.5) +
  geom_density(weibullNormalEagle, mapping = aes(x = weibull), col = "firebrick", lwd = 2.5)

ggplot() + geom_density(eagleData, mapping = aes(x = vertical_steps),
                        col = "dodgerblue4", lwd = 2.5) +
  geom_density(weibullNormalEagle, mapping = aes(x = normal), col = "firebrick", lwd = 2.5)

# And jointly:
ggplot() + #geom_point(data = eagleData, mapping = aes(x = horizontal_steps, y = vertical_steps),
  #alpha = .1) +
  stat_density_2d(data = weibullNormalEagle, mapping = aes(x=weibull, y=normal, fill = ..level..)
                  , geom = "polygon", colour="white") +
  scale_fill_gradient(low = "darkred", high = "firebrick", na.value = NA)


# Idea with clustering
clusters <- stats::kmeans(data.frame(x = eagleData$horizontal_steps,
                                     y = eagleData$vertical_steps), centers = N)$cluster

# eagleDataCluster <- tibble(eagleData, clusters = clusters) %>% mutate(clusters = factor(clusters))
# 
# eagleDataCluster %>% ggplot(aes(x = vertical_steps, fill = clusters)) + geom_density()
# tibble(eagleData, clusters = factor(stats::kmeans(eagleData$vertical_steps, centers = N)$cluster)) %>%
#   ggplot(aes(x = vertical_steps, fill = clusters)) + geom_density()
# 
# eagleDataCluster %>% ggplot(aes(x = horizontal_steps, fill = clusters)) + geom_density()
# tibble(eagleData, clusters = factor(stats::kmeans(eagleData$horizontal_steps, centers = N)$cluster)) %>%
#   ggplot(aes(x = horizontal_steps, fill = clusters)) + geom_density()

clusterK <- numeric(N)
clusterLambda <- numeric(N)
clusterMu <- numeric(N)
clusterSigmaSq <- numeric(N)

for (i in 1:N){
  clusterX <- eagleDataCluster$horizontal_steps[eagleDataCluster$clusters  == i]
  clusterY <- eagleDataCluster$vertical_steps[eagleDataCluster$clusters  == i]
  clusterK[i] <- uniroot(function(k) {sum(clusterX^k * log(clusterX)) / sum(clusterX^k) -
      1/k - mean(log(clusterX))}, c(0.1, 10))$root
  
  clusterLambda[i] <- (mean(clusterX^clusterK[i]))^(1/clusterK[i])
  
  clusterMu[i] <- mean(clusterY)
  clusterSigmaSq[i] <- var(clusterY)
}

eagleClusterFit <- momentuHMM::fitHMM(prepEagleData,
                                      nbStates = N,
                                      dist = list(horizontal_steps = "weibull",
                                                  vertical_steps = "norm"),
                                      Par0 = list(horizontal_steps = c(clusterK, clusterLambda),
                                                  vertical_steps = c(clusterMu, clusterSigmaSq)))


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