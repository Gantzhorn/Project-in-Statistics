Rcpp::sourceCpp("naiveSimRcpp.cpp")
N <- 3
obs1 <- 10000
Gamma1 <- matrix(c(0.5, 0.3, 0.2,
                   0.05, 0.9, 0.05,
                   0.1, 0.3, 0.6), nrow = N, byrow = TRUE)
rowSums(Gamma1) == rep(1, N) # Check if legal matrix
delta1 <- rep(1/N, N)

sum(delta1) == 1

result1 <- fastMarkovChainSampler(obs1, Gamma1, delta1)

# Compute the empirical state probabilities
empirical_probs <- table(result1) / obs1

# Compare with the theoretical stationary distribution
stationary_probs <- eigen(t(Gamma1))$vectors[, 1] / sum(eigen(t(Gamma1))$vectors[, 1])

# Check if the empirical probabilities are close to the stationary distribution
all(abs(empirical_probs - stationary_probs) < 0.05)  # Adjust the tolerance as needed

# Compute the empirical transition matrix
empirical_transitions <- matrix(0, nrow = nrow(Gamma1), ncol = ncol(Gamma1))
for (i in 1:(obs1 - 1)) {
  from_state <- result1[i]
  to_state <- result1[i + 1]
  empirical_transitions[from_state, to_state] <- empirical_transitions[from_state, to_state] + 1
}
empirical_transitions <- empirical_transitions / rowSums(empirical_transitions)

# Compare with the theoretical transition matrix
transition_probs <- Gamma1

# Check if the empirical transition matrix is close to the theoretical matrix
all(abs(empirical_transitions - transition_probs) < 0.05)  # Adjust the tolerance as needed

# Check if the number of occurences of the first value of the chain matches the initial distribution
first_states <- replicate(1000, {
  result <- fastMarkovChainSampler(obs1, Gamma1, delta1)
  result[1]
})

# Count the occurrences of each state as the first entry
occurrences <- table(first_states)

# Compute the empirical initial state distribution
empirical_initial <- occurrences / sum(occurrences)

# Compare with the theoretical initial state distribution
initial_probs <- delta1

# Check if the empirical initial state distribution is close to the theoretical distribution
all(abs(empirical_initial - initial_probs) < 0.05)  # Adjust the tolerance as needed

