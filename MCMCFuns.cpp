
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double nlogDensity (mat data, vec theta, vec mu, mat varInv, int lagsA = 1, int lagsD = 1){
  int T = data.n_rows;
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));

  int lag = max(lagsA, lagsD);
  for(int t = lag; t < T; ++t){
    
    double meanA = 0, meanD = 0;
    for(int i = 1; i <= lagsA; ++i){
      meanA += theta(1 + i) * data(t-i, 0);
    }
    for(int i = 1; i <= lagsD; ++i){
      meanD += theta(1 + lagsA + i) * data(t - i, 1);
    }
 
    double z = pow(data(t, 0) - meanA, 2) / exp(theta(0))  +  pow(data(t, 1) - meanD, 2) / exp(theta(1));
 
    dens += - 0.5 * (theta(0) + theta(1)) - z/2;
    
  }
  return dens;
}







