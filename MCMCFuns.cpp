
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double nlogDensity (mat data, vec theta, vec mu, mat varInv, int lagsA = 1, int lagsD = 1, bool incRho = true){
  int T = data.n_rows;
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  double rho = 0;
  int cons = 1;
  if(incRho){
    rho = 2 / (1 + exp(-theta(2))) - 1;
    cons += 1;
  }

  int lag = max(lagsA, lagsD);
  for(int t = lag; t < T; ++t){
    
    double meanA = 0, meanD = 0;
    for(int i = 1; i <= lagsA; ++i){
      meanA += theta(cons + i) * data(t-i, 0);
    }
    for(int i = 1; i < lagsD; ++i){
      meanD += theta(cons + lagsA + i) * data(t - i, 1);
    }
 
    double z = pow(data(t, 0) - meanA, 2) / exp(theta(0))  +  pow(data(t, 1) - meanD, 2) / exp(theta(1))  -
      2 * rho * (data(t, 0) - meanA) * (data(t, 1) - meanD) / (sqrt(exp(theta(0))) * sqrt(exp(theta(1))));
 
    dens += - 0.5 * (theta(0) + theta(1) + log(1 - rho * rho))  -  z / (2 - 2 * rho * rho);
    
  }
  return dens;
}

double signDouble(double x){
  if(x > 0){
    return 1.0;
  } else if(x < 0) {
    return -1.0;
  } else {
    return 0.0;
  }
}

// [[Rcpp::export]]
double nlogDensityVAR (mat data, vec theta, vec mu, mat varInv, int lagsA = 1, int lagsD = 1, bool incRho = true){
  int T = data.n_rows;
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  double rho = 0;
  int cons = 1;
  if(incRho){
    rho = 2 / (1 + exp(-theta(2))) - 1;
    cons += 1;
  }
  int lag = max(lagsA, lagsD);
  for(int t = lag; t < T; ++t){
    double meanA = 0, meanD = 0;
    for(int i = 1; i <= lagsA; ++i){
      meanA += theta(cons + i) * data(t-i, 0)  +  theta(cons + lagsA + i) * abs(data(t-i, 1)) * signDouble(data(t-i, 0));
    } 
    for(int i = 1; i <= lagsD; ++i){
      meanD += theta(cons + 2 * lagsA + i) * data(t - i, 1)  +  theta(cons + 2 * lagsA + lagsD + i) * abs(data(t-i, 0)) * signDouble(data(t-i, 1));
    }
    
    double z = pow(data(t, 0) - meanA, 2) / exp(theta(0))  +  pow(data(t, 1) - meanD, 2) / exp(theta(1))  -
      2 * rho * (data(t, 0) - meanA) * (data(t, 1) - meanD) / (sqrt(exp(theta(0))) * sqrt(exp(theta(1))));
    
    dens += - 0.5 * (theta(0) + theta(1) + log(1 - rho * rho))  -  z / (2 - 2 * rho * rho);
  }
  return dens;
}







