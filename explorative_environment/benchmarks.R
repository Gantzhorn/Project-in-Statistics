# Run relevant functions in other R files before running this script
slowMarkovChainSampler <- function(obs, Gamma, delta){
  s <- rep(0, obs)
  s[1] <- sample.int(N, size = 1, prob = delta)
  for (i in 2:obs){
    s[i] <- sample.int(N, size = 1, prob = Gamma[s[(i-1)],])
  }
  return(s)
}
plotPath = "/home/anders/Desktop/Skole/Blok 4 - 2023/Project in Statistics/Project-in-Statistics/latex/figures"
# Varying lengths of obs
obs_lengths <- c(500, 1000, 2500, 5000, 10000, 17500, 25000, 35000, 50000, 75000)

# Empty vectors to store the benchmark results
fast_times <- numeric(length(obs_lengths))
fast_times_lower <- numeric(length(obs_lengths))
fast_times_upper <- numeric(length(obs_lengths))
slow_times <- numeric(length(obs_lengths))
slow_times_lower <- numeric(length(obs_lengths))
slow_times_upper <- numeric(length(obs_lengths))

# Benchmark the samplers for each observation length
for (i in seq_along(obs_lengths)) {
  obs_length <- obs_lengths[i]
  
  # Benchmark the fast sampler
  fast <- quantile(microbenchmark(
         fastMarkovChainSampler(obs_length, Gamma, delta),
         times = 50,  # Number of repetitions
         control = list(warmup = 3)  # Warm-up the function before timing
    )$time)[2:4]
  fast_times[i] <- fast[2]
  fast_times_lower[i] <- fast[1]
  fast_times_upper[i] <- fast[3]
  # Benchmark the slow sampler
  slow <- quantile(microbenchmark(
         slowMarkovChainSampler(obs_length, Gamma, delta),
         times = 50,  # Number of repetitions
         control = list(warmup = 3)  # Warm-up the function before timing
    )$time)[2:4]
  slow_times[i] <- slow[2]
  slow_times_lower[i] <- slow[1]
  slow_times_upper[i] <- slow[3]
}

# Create a data frame with the benchmark results
benchmark_df <- data.frame(
  obs_length = obs_lengths,
  sampler = rep(c("C++", "R"), each = length(obs_lengths)),
  time = c(fast_times, slow_times),
  lower_time = c(fast_times_lower, slow_times_lower),
  upper_time = c(fast_times_upper, slow_times_upper)
)
# Plot the benchmark results
bench1 <- ggplot(benchmark_df, aes(x = obs_length, y = time, color = sampler)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) + 
  labs(
    x = "Length of Chain",
    y = "Execution Time (ms)",
    color = "Language"
  ) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = proj_palette) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(1, "cm")
  )
ggsave(filename = "markovChainSamplers.jpeg", plot = bench1, path = plotPath, width = 7, height = 5)
