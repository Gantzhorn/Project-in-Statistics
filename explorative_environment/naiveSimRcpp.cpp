#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix split_and_stack_rcpp(Rcpp::NumericVector x,
                                         int ncol) {
  int len = x.length();
  int nrow = len / ncol;
  Rcpp::NumericVector x_pad = Rcpp::rep(NA_REAL, nrow * ncol);
  x_pad[Rcpp::Range(0, len)] = x;
  Rcpp::NumericMatrix mat(nrow, ncol);
  int k = 0;
  for(int i = 0; i < nrow; i++) {
    for(int j = 0; j < ncol; j++) {
      mat(i, j) = x[k];
      k++;
    }
  }

  return mat;
}

// [[Rcpp::export]]
double jarque_bera_test_Rcpp_Naive(Rcpp::NumericVector x) {
  // Check if there are any NAs in x
  if (is_true(any(is_na(x)))) {
    Rcpp::stop("NAs in x");
  }
  // Compute the test statistic
  int n = x.length();
  double m1 = mean(x);
  double m2 = mean(pow(x - m1, 2));
  double m3 = mean(pow(x - m1, 3));
  double m4 = mean(pow(x - m1, 4));
  double b1 = pow(m3 / pow(m2, 1.5), 2);
  double b2 = m4 / pow(m2, 2);
  double statistic = n * (b1 / 6 + pow(b2 - 3, 2) / 24);
  // Compute the p-value using the chi-square distribution
  double p_value = exp(-statistic/2);
  return p_value;
}

// [[Rcpp::export]]
double jarque_bera_test_Rcpp_optimized(Rcpp::NumericVector x) {
  // Check if there are any NAs in x
  if (is_true(any(is_na(x)))) {
    Rcpp::stop("NAs in x");
  }
  // Compute the test statistic
  int n = x.length();
  double m1 = 0, m2 = 0, m3 = 0, m4 = 0,
    delta1 = 0, delta2 = 0, delta3 = 0, delta4 = 0,
    sNumerator = 0, denominators = 0, kNumerator = 0,
    S_square = 0, K = 0;
  double x_mean = mean(x);
  for (int i = 0; i < n; ++i) {
    delta1 = (x[i]-x_mean);
    delta2 = delta1 * delta1;
    delta3 = delta2 * delta1;
    delta4 = delta2 * delta2;
    sNumerator += delta3;
    denominators += delta2;
    kNumerator += delta4;
  }
  S_square = n*sNumerator*sNumerator/pow(denominators, 3);
  K = n*kNumerator/(denominators*denominators);
  
  double statistic = n * (S_square / 6 + pow(K - 3, 2) / 24);
  // Compute the p-value using the chi-square distribution
  double p_value = exp(-statistic/2);
  return p_value;
}

// [[Rcpp::export]]
Rcpp::IntegerVector fastMarkovChainSampler(int obs,
                                           Rcpp::NumericMatrix Gamma,
                                           Rcpp::NumericVector delta) {
  int nStates = delta.length();
  Rcpp::NumericVector x = Rcpp::runif(obs); // Simulate all comparisons for the chain
  Rcpp::IntegerVector I(obs);
  Rcpp::NumericVector probs(nStates);
  
  // Simulate the first entry in the Markov chain
  double cumulativeSum = 0.0;
  for (int k = 0; k < nStates; k++) {
    cumulativeSum += delta[k];
    probs[k] = cumulativeSum;
  }
  
  for (int j = 0; j < nStates; j++) {
    if (x[0] < probs[j]) {
      I[0] = j + 1;
      break;
    }
  }
  
  // Sample the rest of the states
  for (int i = 1; i < obs; i++) {
    cumulativeSum = 0.0;
    
    // Update weights according to the transition probability matrix
    for (int k = 0; k < nStates; k++) {
      cumulativeSum += Gamma(I[i - 1] - 1, k);
      probs[k] = cumulativeSum;
    }
    
    for (int j = 0; j < nStates; j++) {
      if (x[i] < probs[j]) {
        I[i] = j + 1;
        break;
      }
    }
  }
  return I;
}