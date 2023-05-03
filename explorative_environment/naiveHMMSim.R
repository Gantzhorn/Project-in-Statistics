library(tidyverse)
library(momentuHMM)
theme_set(theme_bw())
proj_palette <- c("#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00",
                         "#CC79A7")

# Define the parameters of the HMM
N <- 3
delta <- c(0.5, 0.25, 0.25)
Gamma <- matrix(c(0.5, 0.25, 0.25,
                  0.25, 0.5, 0.25,
                  0.25, 0.25, 0.5),
            nrow = N, byrow = TRUE) # State transition matrix

# Number of time-steps
T <- 2000
# Parameters for weibull distributions of horizontal distances
k <- c(3, 10, 12)
lambda <- c(2, 4, 6)

# Parameters for normal distributions for vertical distances
mu <- c(-1, 0, 1)
sigma <- c(1, 0.5, 1)

simulate_HMM <- function(T, delta, Gamma, lambda, k, mu, sigma) {
  
  N <- length(delta)
  # Simulate the hidden state sequence
  s <- rep(0, T)
  s[1] <- sample(1:N, size = 1, prob = delta)
  s[2:T] <- sample(1:N, size = T-1, prob = Gamma[s[1:T-1],], replace=TRUE)
  
  # Simulate the observations
  obs_mat <- matrix(nrow = T, ncol = 2)
  
  # X-coordinate
  obs_mat[,1] <- rweibull(T, scale = lambda[s], shape = k[s])
  # Y-coordinate
  obs_mat[,2] <- rnorm(T, mean = mu[s], sd = sigma[s])
  
  # Return the simulated data as a list
  return(list(s = s, obs_mat = obs_mat))
}

simRes <- simulate_HMM(T, delta, Gamma, lambda, k, mu, sigma)

# Return the simulated hidden state sequence and observation sequence
simResult <- tibble(state = factor(simRes$s),
                    x = simRes$obs_mat[,1], y = simRes$obs_mat[,2]) %>%
  tibble::rowid_to_column("time")
# Nice plots
simResult %>% ggplot(aes(x = y, fill = state)) + geom_density()  +
  scale_fill_manual(values = proj_palette)

simResult %>% ggplot(aes(x = x, y = y, color = state)) + geom_point()  +
  scale_color_manual(values = proj_palette)

simResult %>% pivot_longer(cols = c("x", "y"),
                           names_to = "coord_lab", values_to = "coord_val") %>%
  mutate(coord_lab = factor(coord_lab)) %>%
  ggplot(aes(x = coord_val, fill = state)) + geom_density() +
  facet_wrap(~coord_lab, scales = "free") + scale_fill_manual(values = proj_palette)

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

### Accessing elements from the momentuHMM fit object
modDim2$data # data used for fitting
modDim2$mle # Estimate of parameters of observed chain, initial distribution
# & probability transition matrix 

modDim2$CIreal # Estimates and confidence intervals of the above

modDim2$mod # mod estimates, gradient, hessian, numIter, min of neg loglik 
# time to conv. etc.

# play with the objects
# plot confidence intervals for estimates
verticalStepsmat <- cbind(modDim2$CIreal$vertical_steps$est,
                          modDim2$CIreal$vertical_steps$lower,
                          modDim2$CIreal$vertical_steps$upper)

colnames(verticalStepsmat) <- c("estState1", "estState2", "estState3",
                                "lowerState1", "lowerState2", "lowerState3",
                                "upperState1", "upperState2", "upperState3")

verticalSteps <- as_tibble(verticalStepsmat, rownames = "Parameter") %>% 
  pivot_longer(cols = -Parameter, names_to = "Info",values_to = "Value") %>% 
  mutate(State = factor(str_sub(Info, start = -1L, end = -1L)), Parameter = factor(Parameter),
         Type = factor(case_when(str_detect(Info,"estState") ~ "Estimate",
                          str_detect(Info,"lowerState") ~ "Lower",
                          str_detect(Info,"upperState")  ~ "Upper"
                          ))) %>% select(Parameter, Type, State, Value) %>% 
  pivot_wider(names_from = "Type", values_from = "Value") %>% 
  mutate(trueParam = c(mu, sigma))

verticalSteps %>% ggplot(aes(x = State, y = Estimate)) +
  geom_point() + geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
  facet_wrap(~Parameter, scale = "free") +
  geom_point(aes(x = State, y = trueParam), col = "firebrick")
