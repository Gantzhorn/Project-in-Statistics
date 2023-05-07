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
