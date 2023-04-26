library(tidyverse)
library(momentuHMM)
theme_set(theme_bw())
# Define the parameters of the HMM
N <- 3
pi <- c(1/3, 1/3, 1/3)
A <- matrix(c(1/3, 1/3, 1/3,
              1/3, 1/3, 1/3,
              1/3, 1/3, 1/3),
            nrow = N, byrow = TRUE) # State transition matrix

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
  obs_seq[t] <- (hidden_seq[t] == 1)*rbeta(1, shape1 = 0.5, shape2 = 0.5) +
    (hidden_seq[t] == 2)*rnorm(1, mean = 3) + (hidden_seq[t] == 3)*rgamma(1, shape = 5, rate = 2)
}


# Return the simulated hidden state sequence and observation sequence
simResult <- tibble(hidden = factor(hidden_seq), obs = obs_seq)

# plot res
simResult %>% ggplot(aes(x = obs, fill = hidden)) + geom_density()

momentuHMM::fitHMM
