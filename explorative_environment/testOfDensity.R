Rcpp::sourceCpp("/home/anders/Desktop/Skole/Blok 4 - 2023/Project in Statistics/Project-in-Statistics/explorative_environment/manualHMM.cpp")
# Fixed parameters
k1 <- 5
lambda1 <- 10
mu1 <- 1
sigma1 <- 2
covarianceMatrix1 <- matrix(c(1,-0.3,-0.3,1),2,2)
eps1 <- 0.001

sims <- MASS::mvrnorm(1000, mu = c(0,0), Sigma = covarianceMatrix1, empirical = TRUE)

X_prime1 <- qweibull(pnorm(sims[, 1]), shape = k1, scale = lambda1)
Y_prime1 <- sims[,2] * sigma1 + mu1

tibble(x = X_prime1, y = Y_prime1) %>% ggplot(aes(x = x, y = y)) + geom_point()


# Define the grid of points for integration
dx <- 0.01
dy <- 0.01
x <- seq(from = 0.2, to = 20, by = dx)
y <- seq(from = -10, to = 10, by = dy)

# Initialize the sum
sum = 0

# Loop over the grid of points
for(xi in x) {
  for(yi in y) {
    # Evaluate the function at each point
    f_val <- gDensityRcpp(xi, yi, k1, lambda1, mu1, sigma1, covarianceMatrix1, eps1)
    # Add the contribution of this point to the sum
    sum <- sum + f_val * dx * dy
  }
}

# Check if the sum is close to 1
if(abs(sum - 1) < 1e-6) {
  print("The integral of the function over all x and y is close to 1.")
} else {
  print("The integral of the function over all x and y is not close to 1.")
}


xs <- seq(0, 3000, length.out = 100)   # Create a sequence of x values.
ys <- dweibull(xs, shape = 1.218027, scale = 683.82704)  # Compute the density of the Weibull at these x values.

plot(xs, ys, type = "l", 
     main = "Weibull Distribution (shape = 1, scale = 1)", 
     xlab = "x", ylab = "Density")
