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
T <- 300

# Simulate the hidden state sequence
hidden_seq <- rep(0, T)
hidden_seq[1] <- sample(1:N, size = 1, prob = pi)
for (t in 2:T) {
  hidden_seq[t] <- sample(1:N, size = 1, prob = A[hidden_seq[t-1],])
}

# Simulate the observations
obs_mat <- matrix(nrow = T, ncol = 2)
for (t in 1:T) {
  # X-coordinate
  obs_mat[t,1] <- (hidden_seq[t] == 1)*rbeta(1, shape1 = 0.5, shape2 = 0.5) +
    (hidden_seq[t] == 2)*rexp(1, rate = 1) +
    (hidden_seq[t] == 3)*rgamma(1, shape = 5, rate = 4)
  # Y-coordinate
  obs_mat[t,2] <- (hidden_seq[t] == 1)*rnorm(1) +
    (hidden_seq[t] == 2)*rnorm(1, mean = 0, sd = 5) + 
    (hidden_seq[t] == 3)*rnorm(1, mean = 0, sd = 10)
}


# Return the simulated hidden state sequence and observation sequence
simResult <- tibble(state = factor(hidden_seq),
                    x = obs_mat[,1], y = obs_mat[,2],
                    ID = factor(rep(c(1,2,3), each = T/3))) %>% 
  tibble::rowid_to_column("time") %>% 
  mutate(time = case_when(
    ID == "2" ~ time-T/3,
    ID == "3" ~ time-2*(T/3),
    TRUE ~ time
  ))

simResult %>% pivot_longer(cols = c("x", "y"),
                           names_to = "coord_lab", values_to = "coord_val") %>%
  mutate(coord_lab = factor(coord_lab)) %>% 
  ggplot(aes(x = time, y = coord_val)) + geom_line() +
  facet_grid(coord_lab~ID, scales = "free") +
  geom_point(aes(x = time, y = coord_val, color = state))



