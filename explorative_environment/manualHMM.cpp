#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector gDensityRcpp(Rcpp::NumericVector x, Rcpp::NumericVector y,
                                   double k, double lambda, double mu, double sigma, Rcpp::NumericMatrix covariancematrix,
                                   double eps) {
  int n = x.length();
  Rcpp::NumericMatrix covariancematrixInverse(2, 2);
  covariancematrixInverse(0, 0) = covariancematrix(0, 0);
  covariancematrixInverse(1, 1) = covariancematrix(1, 1);
  covariancematrixInverse(0, 1) = -covariancematrix(0, 1);
  covariancematrixInverse(1, 0) = -covariancematrix(1, 0);

  double determinant = covariancematrix(0, 0) * covariancematrix(1, 1) - covariancematrix(0, 1) * covariancematrix(1, 0);
  covariancematrixInverse = covariancematrixInverse / determinant;
  Rcpp::NumericVector exponent = Rcpp::pow(x/lambda, k);
  Rcpp::NumericVector exponential = Rcpp::exp(-exponent);

  Rcpp::NumericVector updateY = 1-exponential;
  for (int i = 0; i < n; i++) {
    if (updateY[i] == 0) {
      updateY[i] += eps;  // Add eps if the entry is equal to 0
    } else if (updateY[i] == 1) {
      updateY[i] -= eps;  // Subtract eps if the entry is equal to 1
    }
  }
  Rcpp::NumericVector gjacobianVec = k * exponent * exponential/(Rcpp::dnorm4(Rcpp::qnorm5(updateY)) * x * sigma + eps);
  Rcpp::NumericVector g1Inverse = (y-mu)/sigma;
  Rcpp::NumericVector g2Inverse = Rcpp::qnorm5(updateY);


  Rcpp::NumericVector quadraticForm = covariancematrixInverse(0,0)*g2Inverse*g2Inverse + 
    covariancematrixInverse(1,1)*g1Inverse*g1Inverse+
    2*covariancematrixInverse(1,0)*g2Inverse*g1Inverse;
  Rcpp::NumericVector density(n);
  
  density = 1.0/(2*3.141592653589793115997963468544185161590576171875*sqrt(determinant))*
    gjacobianVec*Rcpp::exp(-1.0/2.0*quadraticForm);
  return density;
}

/*
 * gDensityVec <- function(x, y, k, lambda, mu, sigma, covarianceMatrix, eps) {
 covMatInverse <- inverseCalculator(covarianceMatrix)
 inputMatrix <- cbind(g1Inverse(y, mu, sigma), g2InverseVec(x, k, lambda, eps))
 result <- 1 / (2 * pi * sqrt(determinantCalculator(covarianceMatrix))) *
 gjacobianVec(x, k, lambda, sigma, eps) *
 exp(-1 / 2 * quadratic_formVec(covMatInverse, inputMatrix))
 return(result)
 }
 */