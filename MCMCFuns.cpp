
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double nlogDensity (mat data, vec theta, vec mu, mat varInv){
  int T = data.n_rows;
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  for(int t = 2; t < T; ++t){
    dens +=  - 0.5 * (theta(0) + theta(1)) - 
      pow(data(t, 0) - theta(2) * data(t-1, 0) - theta(3) * data(t-2, 0), 2) / (2 * exp(theta(0))) -
      pow(data(t, 1) - theta(4) * data(t-1, 1) - theta(5) * data(t-2, 1), 2) / (2 * exp(theta(1)));
  }
  return dens;
}

// [[Rcpp::export]]
double tlogDensity (mat data, vec theta, vec mu, mat varInv){
  int T = data.n_rows;
  double sigSqV = exp(theta(0)), sigSqD = exp(theta(1)), nuV = exp(theta(6)), nuD = exp(theta(7)),
    dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  for(int t = 2; t < T; ++t){
    dens += lgamma((1+nuV)/2)  +  lgamma((1+nuD)/2)  -  lgamma(nuV/2)  -  lgamma(nuD/2) -
      0.5 * (theta(0) + theta(1) + theta(6) + theta(7))  - 
      (nuV+1)/2 * log(1 + pow(data(t, 0) - theta(2) * data(t-1, 0) - theta(3) * data(t-2, 0), 2) / (nuV * sigSqV))  -
      (nuD+1)/2 * log(1 + pow(data(t, 1) - theta(4) * data(t-1, 1) - theta(5) * data(t-2, 1), 2) / (nuD * sigSqD));
  }
  return dens;
}

// [[Rcpp::export]]
double mvtDens (vec theta, vec mu, mat varInv){
  double logDens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  return exp(logDens);
}








