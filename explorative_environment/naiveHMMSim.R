library(tidyverse)
theme_set(theme_bw())
# Define the parameters of the HMM
N <- 3
pi <- c(0.4, 0.3, 0.3)
A <- matrix(c(0.6, 0.2, 0.2, 0.3, 0.5, 0.2, 0.1, 0.4, 0.5), nrow = N, byrow = TRUE) # State transition matrix

# Number of time-steps
T <- 10000

# Simulate the hidden state sequence
hidden_seq <- rep(0, T)
hidden_seq[1] <- sample(1:N, size = 1, prob = pi)
for (t in 2:T) {
  hidden_seq[t] <- sample(1:N, size = 1, prob = A[hidden_seq[t-1],])
}

# Simulate the observations
obs_seq <- rep(0, T)
for (t in 1:T) {
  obs_seq[t] <- (hidden_seq[t] == 1)*runif(1, min = -5, max = -1) +
    (hidden_seq[t] == 2)*rnorm(1, mean = 5) + (hidden_seq[t] == 3)*rgamma(1, shape = 5, rate = 2)
}


# Return the simulated hidden state sequence and observation sequence
simResult <- tibble(hidden = factor(hidden_seq), obs = obs_seq)

# plot res
simResult %>% ggplot(aes(x = obs, fill = hidden)) + geom_density()

