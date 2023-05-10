#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
double jarque_bera_test_Rcpp(Rcpp::NumericVector x) {
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
