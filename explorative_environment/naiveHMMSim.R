library(tidyverse)
library(momentuHMM)
library(Rcpp)
library(profvis)
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
                  0.2, 0.6, 0.2,
                  0.2, 0.2, 0.6),
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
  s <- rep(0, T)
  if(fastSampler){
    s[1] <- sample.int(N, size = 1, prob = delta)
    s[2:T] <- sample.int(N, size = T-1, prob = Gamma[s[1:T-1],], replace=TRUE)
  }
  else{
  s[1] <- sample(1:N, size = 1, prob = delta)
  s[2:T] <- sample(1:N, size = T-1, prob = Gamma[s[1:T-1],], replace=TRUE)
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

simRes <- simulate_HMM(T, delta, Gamma,
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

DecodedStates1 <- viterbi(m = modDim1) #Most likely state-sequence - compare to true state sequence.
mean(DecodedStates1 == simRes$s) # Classification accuracy

mean(DecodedStates1 != simRes$s)

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

plot(modDim2)

DecodedStates2 <- viterbi(m = modDim2) #Most likely state-sequence - compare to true state sequence.


mean(DecodedStates2 == simRes$s) # Classification accuracy

mean(DecodedStates2 != simRes$s)

# Sim from model:
#newSim <- simData(model = modDim2, nbStates = 3, obsPerAnimal = 5000)

### Accessing elements from the momentuHMM fit object
modDim2$data # data used for fitting
modDim2$mle # Estimate of parameters of observed chain, initial distribution
# & probability transition matrix

modDim2$CIreal # Estimates and confidence intervals of the above

modDim2$mod # mod estimates, gradient, hessian, numIter, min of neg loglik
# time to conv. etc.
# 
# # play with the objects
# # plot confidence intervals for estimates
# verticalStepsmat <- cbind(modDim2$CIreal$vertical_steps$est,
#                           modDim2$CIreal$vertical_steps$lower,
#                           modDim2$CIreal$vertical_steps$upper)
# 
# colnames(verticalStepsmat) <- c("estState1", "estState2", "estState3",
#                                 "lowerState1", "lowerState2", "lowerState3",
#                                 "upperState1", "upperState2", "upperState3")
# 
# verticalSteps <- as_tibble(verticalStepsmat, rownames = "Parameter") %>%
#   pivot_longer(cols = -Parameter, names_to = "Info",values_to = "Value") %>%
#   mutate(State = factor(str_sub(Info, start = -1L, end = -1L)), Parameter = factor(Parameter),
#          Type = factor(case_when(str_detect(Info,"estState") ~ "Estimate",
#                           str_detect(Info,"lowerState") ~ "Lower",
#                           str_detect(Info,"upperState")  ~ "Upper"
#                           ))) %>% select(Parameter, Type, State, Value) %>%
#   pivot_wider(names_from = "Type", values_from = "Value") %>%
#   mutate(trueParam = c(mu, sigma))
# 
# verticalSteps %>% ggplot(aes(x = State, y = Estimate)) +
#   geom_point() + geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
#   facet_wrap(~Parameter, scale = "free") +
#   geom_point(aes(x = State, y = trueParam), col = "firebrick")

# All very good let's access fit
pseudoResmodDim2 <- momentuHMM::pseudoRes(modDim2)
pseudoResTibble <- tibble(horizontal_stepsRes = pseudoResmodDim2[[1]],
                          vertical_stepsRes = pseudoResmodDim2[[2]]
                          ) %>% mutate(time = row_number())
# replicate base plots made in this object with ggplot2
momentuHMM::plotPR(modDim2)

# Chain independence:
pseudoResTibble %>% ggplot(aes(x = time, y = horizontal_stepsRes)) +
  geom_line(color = proj_palette[5])

# Histogram of residuals with standard normal distribution on top
pseudoResTibble %>% ggplot(aes(x = vertical_stepsRes)) +
  geom_histogram(aes(y = after_stat(density)),col = "black", fill = proj_palette[6], alpha = 0.8) + 
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = proj_palette[5], linewidth = 1)

# QQ-plot of residuals
pseudoResTibble %>% ggplot(aes(sample = horizontal_stepsRes)) + 
  geom_qq() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = "QQ Plot of Horizontal Steps Pseudo Residuals")

# Make autocorrelation plots
acf_vals <- acf(pseudoResTibble$horizontal_stepsRes, plot = FALSE)
acftib <- tibble(ACF = acf_vals$acf, lag = acf_vals$lag)
ggplot(acftib, aes(x = lag, y = ACF)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(0.025, -0.025), linetype = "dashed", color = proj_palette[3]) +
  geom_segment(aes(xend = lag, yend = 0), color = "black") +
  xlab("Lag") +
  ylab("Autocorrelation")

# More formal test to test fit - jarque bera test on residuals
jarqueBeraTest <- tseries::jarque.bera.test(pseudoResmodDim2$vertical_stepsRes)
jarque_bera_test_Rcpp_Naive(pseudoResmodDim2$vertical_stepsRes)
jarque_bera_test_Rcpp_optimized(pseudoResmodDim2$vertical_stepsRes)
bench_vector <- c(50, 100, 500, 1000, 2500, 10000, 15000, 30000, 50000, 75000, 100000)

# Create a data frame to store the results
results_df <- data.frame(
  length = integer(),
  function1_time = double(),
  function2_time = double()
)

# Loop over the different input sizes and benchmark the functions
for (length in bench_vector) {
  x <- rnorm(length)
  benchmark_results <- microbenchmark(
    tseries::jarque.bera.test(x),
    jarque_bera_test_Rcpp_optimized(x),
    times = 250  # number of times to repeat each function call
  )
  median_times <- aggregate(time/100000000000 ~ expr, benchmark_results, median)
  results_df <- rbind(
    results_df,
    data.frame(
      length = length,
      tseries = median_times$time[1],
      Rcpp = median_times$time[2]
    )
  )
}

# Plot the results
ggplot(results_df, aes(x = length)) +
  geom_line(aes(y = tseries, color = "tseries"), linewidth = .95) +
  geom_line(aes(y = Rcpp, color = "Rcpp"), linewidth = .95) +
  scale_y_log10() + scale_x_log10() + 
  xlab("Input vector length") +
  ylab("Computation time (seconds)") +
  ggtitle("Comparison of tseries- and Rcpp implementation") + 
  scale_color_manual(values = proj_palette)

# Simulate to find bias in fit
M <- 50 # Number of simulations

# Classification
classificationCheck <- numeric(M)

# Bias in estimation
biashorizontal_stepsCheckList <- vector("list", M)
biasvertical_stepsCheckList <- vector("list", M)
biasGammaCheckList <- vector("list", M)

# Confidence intervals
CIhorizontal_stepsCheckList <- vector("list", M)
CIvertical_stepsCheckList <- vector("list", M)
CIGammaCheckList <- vector("list", M)
# Set seed for reproducability
set.seed(1)
for (i in 1:M){
print(paste("Simulating state: ", i))
BiasCheckSim <- simulate_HMM(T, delta, Gamma, dists = c("rweibull", "rnorm"),
                             list(c(k, lambda), c(mu, sigma)))

Prep_dataBiasCheck <- momentuHMM::prepData(data = data.frame(horizontal_steps = BiasCheckSim$obs_mat[,1],
                                           vertical_steps = BiasCheckSim$obs_mat[,2]),
                                           coordNames = NULL
                                           )
modBiasCheck <- fitHMM(data = Prep_dataBiasCheck,
                  nbStates = N,
                  dist = list(horizontal_steps = "weibull",
                              vertical_steps = "norm"),
                  Par0 = list(horizontal_steps = c(k0, lambda0),
                              vertical_steps = c(mu0, sigma0))
)

# Classification
DecodedStatesBiasCheck <- viterbi(m = modBiasCheck) 

classificationCheck[i] <- mean(DecodedStatesBiasCheck == BiasCheckSim$s)

# Bias in estimation
biashorizontal_stepsCheckList[[i]] <- modBiasCheck$mle$horizontal_steps -
  matrix(c(k,lambda), ncol = N, byrow = TRUE)

biasvertical_stepsCheckList[[i]] <- modBiasCheck$mle$vertical_steps -
  matrix(c(mu,sigma), ncol = N, byrow = TRUE)

biasGammaCheckList[[i]] <- modBiasCheck$mle$gamma-Gamma

# Confidence intervals
# Parameters for observed states:
CIGammaCheckList
CIhorizontal_stepsCheckList[[i]] <- (modBiasCheck$CIreal$horizontal_steps$lower <
                                     matrix(c(k,lambda), ncol = N, byrow = TRUE) &
                                     matrix(c(k,lambda), ncol = N, byrow = TRUE) <
                                       modBiasCheck$CIreal$horizontal_steps$upper)

CIvertical_stepsCheckList[[i]] <- (modBiasCheck$CIreal$vertical_steps$lower <
                                   matrix(c(k,lambda), ncol = N, byrow = TRUE) &
                                   matrix(c(mu,sigma), ncol = N, byrow = TRUE) <
                                   modBiasCheck$CIreal$vertical_steps$upper)

CIGammaCheckList[[i]] <- (modBiasCheck$CIreal$gamma$lower < 
                            Gamma & 
                            Gamma <
                            modBiasCheck$CIreal$gamma$upper)
}

# Plot results
# Classification accuracy
tibble(x = classificationCheck) %>% ggplot(aes(x = x)) + 
  geom_density(fill = proj_palette[6], alpha = .85)

# Bias in estimates
# Combinations of row and col are the two last arguments of lapply
tibble(x = unlist(lapply(biashorizontal_stepsCheckList, "[", 1, 1))) %>% ggplot(aes(x = x)) +
  geom_density(fill = proj_palette[6], alpha = .85)

# Bias of transistion probability matrix
tibble(p11 = unlist(lapply(biasGammaCheckList, "[", 1, 1)),
       p12 = unlist(lapply(biasGammaCheckList, "[", 1, 2)),
       p13 = unlist(lapply(biasGammaCheckList, "[", 1, 3)),
       p21 = unlist(lapply(biasGammaCheckList, "[", 2, 1)),
       p22 = unlist(lapply(biasGammaCheckList, "[", 2, 2)),
       p23 = unlist(lapply(biasGammaCheckList, "[", 2, 3)),
       p31 = unlist(lapply(biasGammaCheckList, "[", 3, 1)),
       p32 = unlist(lapply(biasGammaCheckList, "[", 3, 2)),
       p33 = unlist(lapply(biasGammaCheckList, "[", 3, 3))) %>% 
  pivot_longer(cols = everything(), names_to = "Transition", values_to = "Estimate") %>% 
  ggplot(aes(x = Estimate)) + geom_density(fill = proj_palette[5]) + 
  facet_wrap(~Transition, scale = "free")

# Confidence intervals
# k and lambda
sum(unlist(lapply(CIhorizontal_stepsCheckList, "[", 1, 1))/M)
sum(unlist(lapply(CIhorizontal_stepsCheckList, "[", 1, 2))/M)
sum(unlist(lapply(CIhorizontal_stepsCheckList, "[", 1, 3))/M)
sum(unlist(lapply(CIhorizontal_stepsCheckList, "[", 2, 1))/M)
sum(unlist(lapply(CIhorizontal_stepsCheckList, "[", 2, 2))/M)
sum(unlist(lapply(CIhorizontal_stepsCheckList, "[", 2, 3))/M)

# mu and sigma
sum(unlist(lapply(CIvertical_stepsCheckList, "[", 1, 1))/M)
sum(unlist(lapply(CIvertical_stepsCheckList, "[", 1, 2))/M)
sum(unlist(lapply(CIvertical_stepsCheckList, "[", 1, 3))/M)
sum(unlist(lapply(CIvertical_stepsCheckList, "[", 2, 1))/M)
sum(unlist(lapply(CIvertical_stepsCheckList, "[", 2, 2))/M)
sum(unlist(lapply(CIvertical_stepsCheckList, "[", 2, 3))/M)

# Gamma
sum(unlist(lapply(CIGammaCheckList, "[", 1, 1))/M)
sum(unlist(lapply(CIGammaCheckList, "[", 1, 2))/M)
sum(unlist(lapply(CIGammaCheckList, "[", 1, 3))/M)
sum(unlist(lapply(CIGammaCheckList, "[", 2, 1))/M)
sum(unlist(lapply(CIGammaCheckList, "[", 2, 2))/M)
sum(unlist(lapply(CIGammaCheckList, "[", 2, 3))/M)
sum(unlist(lapply(CIGammaCheckList, "[", 3, 1))/M)
sum(unlist(lapply(CIGammaCheckList, "[", 3, 2))/M)
sum(unlist(lapply(CIGammaCheckList, "[", 3, 3))/M)

