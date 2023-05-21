library(tidyverse)
library(momentuHMM)
library(Rcpp)
library(gridExtra)
library(tseries)
library(microbenchmark)
theme_set(theme_bw())


proj_palette <- c("#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00",
                         "#CC79A7")

# Define the parameters of the HMM
N <- 3
delta <- c(0.33, 0.33, 0.33)
Gamma <- matrix(c(0.8, 0.1, 0.1,
                  0.1, 0.8, 0.1,
                  0.1, 0.1, 0.8),
                nrow = N, byrow = TRUE) # State transition matrix

# Number of time-steps
T <- 5000
# Parameters for weibull distributions of horizontal distances
k <- c(3, 10, 12)
lambda <- c(2, 4, 6)

# Parameters for normal distributions for vertical distances
mu <- c(-1, 0, 1)
sigma <- c(1, 0.5, 1)
# Helper function
split_and_stack <- function(x, ncol) {
  # calculate number of rows needed
  nrow <- ceiling(length(x) / ncol)

  # add padding to x to make its length a multiple of ncol
  x_pad <- c(x, rep(NA, nrow * ncol - length(x)))

  # create matrix with desired number of columns and rows
  mat <- matrix(x_pad, nrow = nrow, ncol = ncol, byrow = TRUE)

  return(mat)
}

Rcpp::sourceCpp("naiveSimRcpp.cpp")

simulate_HMM <- function(T = 500,
                         delta,
                         Gamma,
                         dists,
                         params,
                         fastSampler = TRUE
                         ) {

  N <- length(delta)

  # Simulate the hidden state sequence

  if(fastSampler){
    s <- fastMarkovChainSampler(T, Gamma, delta)
  }
  else{
    s <- rep(0, T)
    s[1] <- sample.int(N, size = 1, prob = delta)
    for (i in 2:T){
      s[i] <- sample.int(N, size = 1, prob = Gamma[s[(i-1)],])
    }
  }
  # Simulate the observations
  obs_mat <- matrix(nrow = T, ncol = 2)

  # X-coordinate
  dist1 <- match.fun(dists[1])
  if(fastSampler){params_dist1 <- split_and_stack_rcpp(params[[1]], ncol = N)}
  else {params_dist1 <- split_and_stack(params[[1]], ncol = N)}
  if (nrow(params_dist1) == 1){
    obs_mat[, 1] <- dist1(n = T, params_dist1[1, s])
  }
  else {
    obs_mat[, 1] <- dist1(n = T, params_dist1[1, s], params_dist1[2, s])
  }


  # Y-coordinate
  dist2 <- match.fun(dists[2])
  if(fastSampler){params_dist2 <- split_and_stack_rcpp(params[[2]], ncol = N)}
  else {params_dist2 <- split_and_stack(params[[2]], ncol = N)}
  if (nrow(params_dist2) == 1){
    obs_mat[,2] <- dist2(n = T, params_dist2[1, s])
  }
  else {
    obs_mat[,2] <- dist2(n = T, params_dist2[1, s], params_dist2[2, s])
  }


  # Return the simulated data as a list
  return(list(s = s, obs_mat = obs_mat))
}

simRes <- simulate_HMM(500, delta, Gamma,
                       dists = c("rweibull", "rnorm"),
                       list(c(k, lambda), c(mu, sigma)))
# Fit model
# Init values for mean and variance of gamma and normal
lambda0 <- lambda
k0 <- k
mu0 <- mu
sigma0 <- sigma

# Fit the HMM model with momentuHMM package 1 dimension

# Prepare the data
Prep_data1 <- momentuHMM::prepData(data = data.frame(vertical_steps = simRes$obs_mat[,2]),
                                  coordNames = NULL)

modDim1 <- fitHMM(data = Prep_data1,
              nbStates = N,
              dist = list(vertical_steps = "norm"),
              Par0 = list(vertical_steps = c(mu0, sigma0))
)
plot(modDim1)

classificationAccuracy <- function(decodedstates, toyData){
  return(mean(decodedstates == toyData$s))
}

DecodedStates1 <- viterbi(m = modDim1) #Most likely state-sequence - compare to true state sequence.
classificationAccuracy(DecodedStates1, simRes) # Classification accuracy

# Fit the HMM model with the momentuHMM package 2 dimensions

# Prepare the data
Prep_data2 <- momentuHMM::prepData(data = data.frame(horizontal_steps = simRes$obs_mat[,1],
                                  vertical_steps = simRes$obs_mat[,2]),
                                  coordNames = NULL)

modDim2 <- fitHMM(data = Prep_data2,
              nbStates = N,
              dist = list(horizontal_steps = "weibull",
                          vertical_steps = "norm"),
              Par0 = list(horizontal_steps = c(k0, lambda0),
                          vertical_steps = c(mu0, sigma0))
              )

DecodedStates2 <- viterbi(m = modDim2) #Most likely state-sequence - compare to true state sequence.
classificationAccuracy(DecodedStates2, simRes) # Classification accuracy

# All very good let's access fit
pseudoResmodDim2 <- momentuHMM::pseudoRes(modDim2)

pseudoResAssessmentHMM <- function(pseudoResidual){
# Number of columns in grid
  if(length(pseudoResidual) > 1){
    gridColsNum <- 2
  }
  else{gridColsNum <- 1}
# Get tibble with what we need to do plots
pseudoResTibble <- as_tibble(do.call(cbind, pseudoResidual)) %>% mutate(time = row_number())

# Chain independence:
chainIndependence  <- purrr::map(names(pseudoResTibble)[1:(length(names(pseudoResTibble))-1)],
  function(column_name) {ggplot(pseudoResTibble, aes(x = time, y = .data[[column_name]])) +
    geom_line(color = proj_palette[5]) +
    labs(title = column_name)
})

chainIndependencePlot <- do.call(grid.arrange, c(chainIndependence, ncol = gridColsNum))

# Histogram of residuals with standard normal distribution on top
histogramNormality <- purrr::map(names(pseudoResTibble)[1:(length(names(pseudoResTibble))-1)],
           function(column_name) {
             ggplot(pseudoResTibble, aes(x = .data[[column_name]])) +
               geom_histogram(aes(y = after_stat(density)),col = "black", fill = proj_palette[6], alpha = 0.8) + 
               stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = proj_palette[5], linewidth = 1)
})
chainIndependencePlot <- do.call(grid.arrange, c(histogramNormality, ncol = gridColsNum))

# QQ-plot of residuals
qqNormality <- purrr::map(names(pseudoResTibble)[1:(length(names(pseudoResTibble))-1)],
           function(column_name) {
             ggplot(pseudoResTibble, aes(sample = .data[[column_name]])) +
               geom_qq() +
               geom_abline(intercept = 0, slope = 1, color = proj_palette[3]) + xlab(as.character(column_name))
           })

qqNormalityPlot <- do.call(grid.arrange, c(qqNormality, ncol = gridColsNum))

# Make autocorrelation plots
acf_plots <- map(names(pseudoResTibble)[1:(length(names(pseudoResTibble))-1)], function(column_name) {
  acf_result <- acf(pseudoResTibble[[column_name]], plot = FALSE)
  
  # Create a data frame of lag and ACF values
  acf_data <- data.frame(lag = acf_result$lag, ACF = acf_result$acf)
  
  # Plot the ACF with lag
  ggplot(acf_data, aes(x = lag, y = ACF)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(0.025, -0.025), linetype = "dashed", color = proj_palette[3], lwd = 1) +
    geom_segment(aes(xend = lag, yend = 0), color = "black") +
    xlab("Lag") +
    ylab("Autocorrelation") +
    labs(title = column_name)
})

# Combine the ACF plots using grid.arrange
acfGridPlot <- gridExtra::grid.arrange(grobs = acf_plots, ncol = gridColsNum)

return(list(pseudoResiduals = pseudoResTibble,
            chainIndependence = chainIndependencePlot,
            qqNormality = qqNormalityPlot,
            acfGrid = acfGridPlot))
}
# More formal test to test fit - jarque bera test on residuals
jarqueBeraTests <- function(pseudoResidual){
  jarqueBeraVector <- numeric(length(pseudoResidual))
  for (i in seq_along(jarqueBeraVector)) {
    jarqueBeraVector[i] <- jarque_bera_test_Rcpp_optimized(pseudoResidual[[i]])
  }
  return(jarqueBeraVector)
}

simStudyHMM <- function(numSim,
                        observationsPerSim,
                        delta,
                        Gamma,
                        simDistributions,
                        simdistributionArguments,
                        fitDistributions,
                        fitDistributionArguments,
                        seed = FALSE){
# Make sure we know the number of states
N <- length(delta)
# Initiliaze vectors and lists to hold data
# Classification
classificationCheck <- numeric(numSim)

# Bias in estimation
biashorizontal_stepsCheckList <- vector("list", numSim)
biasvertical_stepsCheckList <- vector("list", numSim)
biasGammaCheckList <- vector("list", numSim)

# Confidence intervals
CIhorizontal_stepsCheckList <- vector("list", numSim)
CIvertical_stepsCheckList <- vector("list", numSim)
CIGammaCheckList <- vector("list", numSim)
# Set seed  for reproducability if chosen to
if(seed){set.seed(1)}


for (i in 1:numSim){
print(paste("Simulating: ", i, " dataset"))
BiasCheckSim <- simulate_HMM(observationsPerSim, delta, Gamma,
                             simDistributions, simdistributionArguments)
# We will assume that we are working with a two state model.
Prep_dataBiasCheck <- momentuHMM::prepData(data = data.frame(horizontal_steps =
                                           BiasCheckSim$obs_mat[,1],
                                           vertical_steps = BiasCheckSim$obs_mat[,2]),
                                           coordNames = NULL
                                           )
print(paste("Fitting HMM number: ", i))
modBiasCheck <- suppressMessages(
  fitHMM(data = Prep_dataBiasCheck,
                       nbStates = N,
                       dist = list(horizontal_steps = fitDistributions[1],
                                   vertical_steps = fitDistributions[2]),
                       Par0 = list(horizontal_steps = fitDistributionArguments[[1]],
                                   vertical_steps = fitDistributionArguments[[2]]))
)
print(paste("Fitting done: ", i))

# Classification accuracy:
classificationCheck[i] <- classificationAccuracy(decodedstates =  viterbi(m = modBiasCheck),
                                                 toyData = BiasCheckSim)
# Bias in estimation
biashorizontal_stepsCheckList[[i]] <- modBiasCheck$mle$horizontal_steps -
  matrix(simdistributionArguments[[1]], ncol = N, byrow = TRUE)

biasvertical_stepsCheckList[[i]] <- modBiasCheck$mle$vertical_steps -
  matrix(simdistributionArguments[[2]], ncol = N, byrow = TRUE)

biasGammaCheckList[[i]] <- modBiasCheck$mle$gamma-Gamma

# Confidence intervals
# Parameters for observed states:
CIhorizontal_stepsCheckList[[i]] <- (modBiasCheck$CIreal$horizontal_steps$lower <
                                       matrix(simdistributionArguments[[1]], ncol = N,
                                              byrow = TRUE) &
                                       modBiasCheck$CIreal$horizontal_steps$upper >
                                       matrix(simdistributionArguments[[1]], ncol = N,
                                              byrow = TRUE))

CIvertical_stepsCheckList[[i]] <- (modBiasCheck$CIreal$vertical_steps$lower <
                                     matrix(simdistributionArguments[[2]], ncol = N,
                                            byrow = TRUE) &
                                     modBiasCheck$CIreal$vertical_steps$upper > 
                                     matrix(simdistributionArguments[[2]], ncol = N,
                                            byrow = TRUE))

CIGammaCheckList[[i]] <- (modBiasCheck$CIreal$gamma$lower < 
                            Gamma & 
                            Gamma <
                            modBiasCheck$CIreal$gamma$upper)
}
return(
  list(numSim = numSim,
       classificationCheck = classificationCheck,
       biashorizontal_stepsCheckList = biashorizontal_stepsCheckList,
       biasvertical_stepsCheckList = biasvertical_stepsCheckList,
       biasGammaCheckList = biasGammaCheckList,
       CIhorizontal_stepsCheckList = CIhorizontal_stepsCheckList,
       CIvertical_stepsCheckList = CIvertical_stepsCheckList,
       CIGammaCheckList = CIGammaCheckList
       )
  )
}

dummy <- simStudyHMM(numSim = 10, observationsPerSim =  500, delta = delta, Gamma = Gamma, 
            simDistributions = c("rweibull", "rnorm"),
            simdistributionArguments = list(c(k, lambda), c(mu, sigma)),
            fitDistributions = c("weibull", "norm"),
            fitDistributionArguments = list(c(k0, lambda0), c(mu0, sigma0)), seed = FALSE)
# Plot results
# Classification accuracy
tibble(x = dummy$classificationCheck) %>% ggplot(aes(x = x)) + 
  geom_density(fill = proj_palette[6], alpha = .85)

# Bias in estimates
# Combinations of row and col are the two last arguments of lapply
tibble(x = unlist(lapply(dummy$biashorizontal_stepsCheckList, "[", 1, 1))) %>% ggplot(aes(x = x)) +
  geom_density(fill = proj_palette[6], alpha = .85)

# Bias of transistion probability matrix
tibble(p11 = unlist(lapply(dummy$biasGammaCheckList, "[", 1, 1)),
       p12 = unlist(lapply(dummy$biasGammaCheckList, "[", 1, 2)),
       p13 = unlist(lapply(dummy$biasGammaCheckList, "[", 1, 3)),
       p21 = unlist(lapply(dummy$biasGammaCheckList, "[", 2, 1)),
       p22 = unlist(lapply(dummy$biasGammaCheckList, "[", 2, 2)),
       p23 = unlist(lapply(dummy$biasGammaCheckList, "[", 2, 3)),
       p31 = unlist(lapply(dummy$biasGammaCheckList, "[", 3, 1)),
       p32 = unlist(lapply(dummy$biasGammaCheckList, "[", 3, 2)),
       p33 = unlist(lapply(dummy$biasGammaCheckList, "[", 3, 3))) %>% 
  pivot_longer(cols = everything(), names_to = "Transition", values_to = "Estimate") %>% 
  ggplot(aes(x = Estimate)) + geom_density(fill = proj_palette[5]) + 
  facet_wrap(~Transition, scale = "free")

# Confidence intervals
# k and lambda
sum(unlist(lapply(dummy$CIhorizontal_stepsCheckList, "[", 1, 1))/dummy$numSim)
sum(unlist(lapply(dummy$CIhorizontal_stepsCheckList, "[", 1, 2))/dummy$numSim)
sum(unlist(lapply(dummy$CIhorizontal_stepsCheckList, "[", 1, 3))/dummy$numSim)
sum(unlist(lapply(dummy$CIhorizontal_stepsCheckList, "[", 2, 1))/dummy$numSim)
sum(unlist(lapply(dummy$CIhorizontal_stepsCheckList, "[", 2, 2))/dummy$numSim)
sum(unlist(lapply(dummy$CIhorizontal_stepsCheckList, "[", 2, 3))/dummy$numSim)

# mu and sigma
sum(unlist(lapply(dummy$CIvertical_stepsCheckList, "[", 1, 1))/dummy$numSim)
sum(unlist(lapply(dummy$CIvertical_stepsCheckList, "[", 1, 2))/dummy$numSim)
sum(unlist(lapply(dummy$CIvertical_stepsCheckList, "[", 1, 3))/dummy$numSim)
sum(unlist(lapply(dummy$CIvertical_stepsCheckList, "[", 2, 1))/dummy$numSim)
sum(unlist(lapply(dummy$CIvertical_stepsCheckList, "[", 2, 2))/dummy$numSim)
sum(unlist(lapply(dummy$CIvertical_stepsCheckList, "[", 2, 3))/dummy$numSim)

# Gamma
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 1, 1))/dummy$numSim)
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 1, 2))/dummy$numSim)
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 1, 3))/dummy$numSim)
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 2, 1))/dummy$numSim)
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 2, 2))/dummy$numSim)
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 2, 3))/dummy$numSim)
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 3, 1))/dummy$numSim)
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 3, 2))/dummy$numSim)
sum(unlist(lapply(dummy$CIGammaCheckList, "[", 3, 3))/dummy$numSim)

