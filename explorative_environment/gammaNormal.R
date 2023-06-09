library(tidyverse)
library(ggExtra)
library(momentuHMM)
library(gridExtra)
covMat <- matrix(c(1, 0.7,
                   0.7, 1), ncol = 2,
                 byrow = TRUE)


bivariateNorm = MASS::mvrnorm(n=50000, mu=c(0, 0),
                              Sigma=covMat, empirical=TRUE)
X = bivariateNorm[, 1]
Y = bivariateNorm[, 2]
cor(X,Y)

gammaNormal <- tibble(gamma = qgamma(pnorm(bivariateNorm[,2], 0, sqrt(covMat[2,2])),
                                     shape = 5, rate = 0.5),
                      normal = bivariateNorm[,1]*2)

# Recall E(X) = alpha/beta and Var(X) = alpha/beta^2
beta_fn <- function(beta){return(mean(gammaNormal$gamma)/beta-var(gammaNormal$gamma))}
beta_hat <- uniroot(beta_fn, interval = c(0.1, 30))$root
alpha_hat <- mean(gammaNormal$gamma)*beta_hat

ggplot(tibble(x = 1:30),aes(x))+
  stat_function(fun=beta_fn) + geom_point(data = tibble(x = beta_hat, y = beta_fn(beta_hat)),
                                          aes(x = x, y = y))

ggplot() +
  geom_density(data = gammaNormal, aes(x = gamma), 
               col = "black", fill = "dodgerblue1") +
  stat_function(fun = dgamma, args = list(shape = alpha_hat, rate = beta_hat), color = "red", size = 1)

ggplot() +
  geom_density(data = gammaNormal, aes(x = normal), 
               col = "black", fill = "dodgerblue1") +
  stat_function(fun = dnorm,
                args = list(mean= mean(gammaNormal$normal), sd = sd(gammaNormal$normal)),
                color = "red", size = 1)

p <- gammaNormal %>% ggplot() + geom_point(aes(x = normal, y = gamma), col = "firebrick", size = 0.5) +
  xlab("Normal") + ylab("Gamma")

ggMarginal(p, fill = "dodgerblue4")


# Actual data fitting
GoldenEagle_FinalData_Dryad_FE <- read_csv(
  "~/Desktop/Skole/Blok 4 - 2023/Project in Statistics/Project-in-Statistics/rmarkdown/GoldenEagle_FinalData_Dryad_FE.csv")

eagleData <- GoldenEagle_FinalData_Dryad_FE %>% group_by(Segment_ID) %>%
  mutate(vertical_steps = Altitude - lag(Altitude, default = first(Altitude))) %>% 
  rename(horizontal_steps = Step_length) %>% 
  slice(2:n())

# The idea:
k_hat <- uniroot(function(k) {sum(eagleData$horizontal_steps^k * log(eagleData$horizontal_steps)) / sum(eagleData$horizontal_steps^k) - 1/k - mean(log(eagleData$horizontal_steps))},
                    c(0.1, 10))$root
lambda_hat <- (mean(eagleData$horizontal_steps^k_hat))^(1/k_hat)

mu_hat <- mean(eagleData$vertical_steps)
sigmasq_hat <- var(eagleData$vertical_steps)

covMatEagle <- matrix(c(1, -0.3,
                   -0.3, 1), ncol = 2,
                 byrow = TRUE)


bivariateNormEagle = MASS::mvrnorm(n=50000, mu=c(0, 0),
                              Sigma=covMatEagle, empirical=TRUE)
X = bivariateNormEagle[, 1]
Y = bivariateNormEagle[, 2]

weibullNormalEagle <- tibble(weibull = qweibull(pnorm(bivariateNormEagle[,2], 0, sqrt(covMatEagle[2,2])),
                                     shape = k_hat, scale = lambda_hat),
                      normal = bivariateNormEagle[,1]*sqrt(sigmasq_hat)+mu_hat)
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
localDecoding <- momentuHMM::stateProbs(eagleFit1)

max_indices <- max.col(a)

statesForSim <- rep(0, nrow(eagleFit1$data))

statesForSim[1] <- sample.int(N, 1, prob = eagleFit1$mle$delta[1,])

for(i in 2:nrow(eagleFit1$data)){
  statesForSim[i] <- sample.int(N, 1, prob = eagleFit1$mle$gamma[statesForSim[(i-1)], ])
}
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

eagleDataCluster <- tibble(eagleData, clusters = clusters) %>% mutate(clusters = factor(clusters))
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
