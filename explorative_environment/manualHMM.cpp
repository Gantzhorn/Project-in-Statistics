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

// [[Rcpp::export]]
double negative_log_likelihood_bivariate_weibull_normal_Rcpp(Rcpp::NumericVector step1,
                                                                Rcpp::NumericVector step2,
                                                                Rcpp::NumericVector mu,
                                                                Rcpp::NumericVector sigma,
                                                                Rcpp::NumericVector shape,
                                                                Rcpp::NumericVector scale,
                                                                Rcpp::NumericMatrix Gamma,
                                                                Rcpp::NumericMatrix covarianceMatrix,
                                                                Rcpp::NumericVector delta,
                                                                double eps = 0.001){
  int T = step1.length();
  int N = delta.length();
  Rcpp::NumericMatrix all_probs(T,3);
  for (int i = 0; i < N; i++) {
      all_probs(Rcpp::_ , i) = gDensityRcpp(step1, step2, shape[i], scale[i], mu[i], sigma[i], covarianceMatrix, eps);
  }
  
  Rcpp::NumericVector v(N);
  double llk = 0.0;
  
  // Initialization
  for (int i = 0; i < N; i++) {
    v[i] = delta[i] * all_probs(0, i);
  }
  
  // Loop over time steps
  for (int t = 1; t < T; t++) {
    // Matrix-vector multiplication
    for (int i = 0; i < N; i++) {
      double sumValue = 0.0;
      for (int j = 0; j < N; j++) {
        sumValue += v[j] * Gamma(j, i);
      }
      v[i] = sumValue * all_probs(t, i);
    }
    
    // Log-sum scaling
    double sumV = sum(v);
    llk += log(sumV);
    v = v / sumV;
  }
  return -llk;
}


