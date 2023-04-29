library(tidyverse)
library(momentuHMM)
theme_set(theme_bw())
# Define the parameters of the HMM
N <- 3
delta <- c(0.5, 0.25, 0.25)
Gamma <- matrix(c(0.8, 0.1, 0.1,
              0.1, 0.8, 0.1,
              0.05, 0.05, 0.9),
            nrow = N, byrow = TRUE) # State transition matrix

# Number of time-steps
T <- 500

# Simulate the hidden state sequence
s <- rep(0, T)
s[1] <- sample(1:N, size = 1, prob = delta)
s[2:T] <- sample(1:N, size = T-1, prob = Gamma[s[1:T-1],], replace=TRUE)
# Simulate the observations
obs_mat <- matrix(nrow = T, ncol = 2)

# Parameters for gamma distributions of horizontal distances
alpha <- c(1, 2, 2)
beta <- c(2, 1, 0.25)

# Parameters for normal distributions for vertical distances
mu <- c(-1, 0, 1)
sigma <- c(0.5, 0.25, 0.5)

# X-coordinate
obs_mat[,1] <- rgamma(T, shape = alpha[s], rate = beta[s])
# Y-coordinate
obs_mat[,2] <- rnorm(T, mean = mu[s], sd = sigma[s])



# Return the simulated hidden state sequence and observation sequence
# simResult <- tibble(state = factor(s),
#                     x = obs_mat[,1], y = obs_mat[,2],
#                     ID = factor(rep(c(1,2,3), each = T/3))) %>% 
#   tibble::rowid_to_column("time") %>% 
#   mutate(time = case_when(
#     ID == "2" ~ time-T/3,
#     ID == "3" ~ time-2*(T/3),
#     TRUE ~ time
#   ))
# Nice plots
# simResult %>% pivot_longer(cols = c("x", "y"),
#                            names_to = "coord_lab", values_to = "coord_val") %>%
#   mutate(coord_lab = factor(coord_lab)) %>% 
#   ggplot(aes(x = time, y = coord_val)) + geom_line() +
#   facet_grid(coord_lab~ID, scales = "free") +
#   geom_point(aes(x = time, y = coord_val, color = state))

# Fit model
# Init values for mean and variance of gamma and normal
alpha0 <- c(1, 2, 2)
beta0 <- c(2, 1, 0.25)

mu0 <- c(-1, 0, 1)
sigma0 <- c(0.5, 0.25, 0.5)

# Fit the HMM model with momentuHMM package 1 dimension

# Prepare the data
Prep_data <- momentuHMM::prepData(data = data.frame(vertical_steps = obs_mat[,2]),
                                  coordNames = NULL)

mod <- fitHMM(data = Prep_data,
              nbStates = N,
              dist = list(vertical_steps = "norm"),
              Par0 = list(vertical_steps = c(mu0, sigma0))
)
plot(mod)


DecodedStates <- viterbi(m = mod) #Most likely state-sequence - compare to true state sequence.
mean(DecodedStates == s) # Classification accuracy

mean(DecodedStates != s)

# Fit the HMM model with the momentuHMM package 2 dimensions

# Prepare the data
Prep_data <- momentuHMM::prepData(data = data.frame(horizontal_steps = obs_mat[,1],
                                  vertical_steps = obs_mat[,2]),
                                  coordNames = NULL)

# Prep_data <- momentuHMM::prepData(data = data.frame(vertical_steps = obs_mat[,2]),
#                                   coordNames = NULL)

mod <- fitHMM(data = Prep_data,
              nbStates = N,
              dist = list(horizontal_steps = "gamma",
                          vertical_steps = "norm"),
              Par0 = list(horizontal_steps = c(alpha0, beta0),
                          vertical_steps = c(mu0, sigma0))
              )

plot(mod)

DecodedStates <- viterbi(m = mod) #Most likely state-sequence - compare to true state sequence.


mean(DecodedStates == s) # Classification accuracy

mean(DecodedStates != s)
