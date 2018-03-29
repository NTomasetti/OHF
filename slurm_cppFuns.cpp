// [[Rcpp::depends(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <boost/math/distributions.hpp> 

using namespace arma;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

// Sobol generation and shuffling

// [[Rcpp::export]]
mat sobol_points(int N, int D) {
  std::ifstream infile("new-joe-kuo-6.21201" ,ios::in);
  if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
    exit(1);
  }
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  
  // L = max number of bits needed 
  unsigned L = (unsigned)ceil(log((double)N)/log(2.0)); 
  
  // C[i] = index from the right of the first zero bit of i
  unsigned *C = new unsigned [N];
  C[0] = 1;
  for (unsigned i=1;i<=N-1;i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // POINTS[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  mat POINTS(N, D, fill::zeros);
  
  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  unsigned *V = new unsigned [L+1]; 
  for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1
  
  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  unsigned *X = new unsigned [N];
  X[0] = 0;
  for (unsigned i=1;i<=N-1;i++) {
    X[i] = X[i-1] ^ V[C[i-1]];
    POINTS(i, 0) = (double)X[i]/pow(2.0,32); // *** the actual points
    //        ^ 0 for first dimension
  }
  
  // Clean up
  delete [] V;
  delete [] X;
  
  
  // ----- Compute the remaining dimensions -----
  for (unsigned j=1;j<=D-1;j++) {
    
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i=1;i<=s;i++) infile >> m[i];
    
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    unsigned *V = new unsigned [L+1];
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
      for (unsigned i=s+1;i<=L;i++) {
        V[i] = V[i-s] ^ (V[i-s] >> s); 
        for (unsigned k=1;k<=s-1;k++) 
          V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
      }
    }
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [N];
    X[0] = 0;
    for (unsigned i=1;i<=N-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
      POINTS(i, j) = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
    }
    
    // Clean up
    delete [] m;
    delete [] V;
    delete [] X;
  }
  delete [] C;
  
  return POINTS;
}

// [[Rcpp::export]]
mat shuffle(mat sobol){
  int N = sobol.n_rows;
  int D = sobol.n_cols;
  mat output(N, D, fill::zeros);
  // draw a random rule of: switch 1 and 0  /  do not switch for each binary digit.
  vec rule = randu<vec>(16);
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < D; ++j){
      // grab element of the sobol sequence
      double x = sobol(i, j);
      // convert to a binary representation
      uvec binary(16, fill::zeros);
      for(int k = 1; k < 17; ++k){
        if(x > pow(2, -k)){
          binary(k-1) = 1;
          x -= pow(2, -k);
        }
      }
      // apply the transform of tilde(x_k) = x_k + a_k mod 2, where a_k = 1 if rule_k > 0.5, 0 otherwise
      for(int k = 0; k < 16; ++k){
        if(rule(k) > 0.5){
          binary(k) = (binary(k) + 1) % 2;
        }
      }
      // reconstruct base 10 number from binary representation
      for(int k = 0; k < 16; ++k){
        if(binary(k) == 1){
          output(i, j) += pow(2, -(k+1));
        }
      }
    }
  }
  return output;
}

// Various VB gradient estimators

// Single Component Prior / Single Component Approximation

struct SPSA {
  const mat data;
  const vec epsilon;
  const vec mean;
  const mat Linv;
  const int lagsA;
  const int lagsD;
  SPSA(const mat& dataIn, const vec& epsIn, const vec& meanIn, const mat& LinvIn, const int& lagsAIn, const int& lagsDIn) :
    data(dataIn), epsilon(epsIn), mean(meanIn), Linv(LinvIn), lagsA(lagsAIn), lagsD(lagsDIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs; using std::max;
    int N = data.n_rows;
    int dim = 2 + lagsA + lagsD;
    // Create theta as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> theta(dim);
    for(int i = 0; i < dim; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(dim*(i+1) + j) * epsilon(j);
      }
    }
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    // Evaluate log(p(theta))
    T prior = 0;
    Matrix<T, Dynamic, 1> kernel(dim);
    kernel.fill(0);
    for(int i = 0; i < dim; ++i){
      for(int j = 0; j <= i; ++j){
        kernel(i) += (theta(j) - mean(j)) * Linv(i, j);
      }
      prior += - 0.5 * pow(kernel(i), 2);
    }
    
    // Evaluate Log Det J
    T logdetJ = theta(0)  +  theta(1);
    for(int i = 0; i < dim; ++i){
      logdetJ += log(fabs(lambda((dim+1)*i+dim)));
    }
    // Evaluate log likelihood
    T logLik = 0;
    Matrix<T, Dynamic, Dynamic> loglikKernel(N, 2);
    int lags = max(lagsA, lagsD);
    for(int t = lags; t < N; ++t){
      loglikKernel(t, 0) = data(t, 0);
      loglikKernel(t, 1) = data(t, 1);
      for(int i = 1; i <= lagsA; ++i){
        loglikKernel(t, 0) -= data(t-i, 0) * theta(1 + i);
      }
      for(int i = 1; i <= lagsD; ++i){
        loglikKernel(t, 1) -= data(t-i, 1) * theta(1 + lagsA + i);
      }
      
      logLik += - 0.5 * theta(0) - 0.5 * theta(1) - 
        pow(loglikKernel(t, 0), 2) / (2 * sigSqV) -
        pow(loglikKernel(t, 1), 2) / (2 * sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List singlePriorSingleApprox(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec mean, mat Linv, int lagsA = 1, int lagsD = 1){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1>  grad(dim);
  // Autodiff
  SPSA p(data, epsilon, mean, Linv, lagsA, lagsD);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

// Mixture Prior / Single Approx

struct MPSA {
  const mat data;
  const vec epsilon;
  const mat mean;
  const cube Linv;
  const vec dets;
  const vec weights;
  const int lagsA;
  const int lagsD;
  MPSA(const mat& dataIn, const vec& epsIn, const mat& meanIn, const cube& LinvIn, const vec& detsIn, const vec& weightsIn, const int& lagsAIn, const int& lagsDIn) :
    data(dataIn), epsilon(epsIn), mean(meanIn), Linv(LinvIn), dets(detsIn), weights(weightsIn), lagsA(lagsAIn), lagsD(lagsDIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs; using std::max;
    int N = data.n_rows;
    int dim = 2 + lagsA + lagsD;
    int mix = weights.n_rows;
    // Create theta as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> theta(dim);
    for(int i = 0; i < dim; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(dim*(i+1) + j) * epsilon(j);
      }
    }
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    
    // Evaluate log(p(theta))
    T prior = 0;
    Matrix<T, Dynamic, Dynamic> kernel(dim, mix);
    kernel.fill(0);
    Matrix<T, Dynamic, 1> exponents(mix);
    exponents.fill(0);
    for(int k = 0; k < mix; ++k){
      for(int i = 0; i < dim; ++i){
        for(int j = 0; j <= i; ++j){
          kernel(i, k) += (theta(j) - mean(j, k)) * Linv(i, j, k);
        }
        exponents(k) += -0.5 * pow(kernel(i, k), 2);
      }
      prior += weights(k) * pow(6.283185, -3) * dets(k) * exp(exponents(k));
    }
    prior = log(prior);
    
    // Evaluate Log Det J
    T logdetJ = theta(0)  +  theta(1);
    for(int i = 0; i < dim; ++i){
      logdetJ += log(fabs(lambda((dim+1)*i+dim)));
    }
    // Evaluate log likelihood
    T logLik = 0;
    Matrix<T, Dynamic, Dynamic> loglikKernel(N, 2);
    int lags = max(lagsA, lagsD);
    for(int t = lags; t < N; ++t){
      loglikKernel(t, 0) = data(t, 0);
      loglikKernel(t, 1) = data(t, 1);
      for(int i = 1; i <= lagsA; ++i){
        loglikKernel(t, 0) -= data(t-i, 0) * theta(1 + i);
      }
      for(int i = 1; i <= lagsD; ++i){
        loglikKernel(t, 1) -= data(t-i, 1) * theta(1 + lagsA + i);
      }
      
      logLik += - 0.5 * theta(0) - 0.5 * theta(1) - 
        pow(loglikKernel(t, 0), 2) / (2 * sigSqV) -
        pow(loglikKernel(t, 1), 2) / (2 * sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List mixPriorSingleApprox(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, mat mean, cube Linv, vec dets, vec weights, int lagsA = 1, int lagsD = 1){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1>  grad(dim);
  // Autodiff
  MPSA p(data, epsilon, mean, Linv, dets, weights, lagsA, lagsD);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

// Pieces for the mixture approximations are shared between both prior specifications

struct mixQlogdens {
  const vec theta;
  const int mix;
  mixQlogdens(const vec& thetaIn, const int& mixIn) :
    theta(thetaIn), mix(mixIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt;
    
    int dim = theta.n_rows;

    Matrix<T, Dynamic, 1> dets(mix);
    for(int k = 0; k < mix; ++k){
      dets(k) = exp(lambda(dim*k + dim*mix));
      for(int i = 1; i < dim; ++i){
        dets(k) *= exp(lambda(dim*mix + dim*k + i));
      }
      dets(k) = 1.0 / dets(k);
    }
    
    Matrix<T, Dynamic, 1> kernel(mix);
    kernel.fill(0);
    
    for(int k = 0; k < mix; ++k){
      for(int i = 0; i < dim; ++i){
        kernel(k) += pow((theta(i) - lambda(k*dim + i)) / exp(lambda(dim*mix + dim*k + i)), 2);
      }
    }
    
    Matrix<T, Dynamic, 1> pi(mix);
    T sumExpZ = 0;
    for(int k = 0; k < mix; ++k){
      pi(k) = exp(lambda(2*dim*mix + k));
      sumExpZ += pi(k);
    }
    pi /= sumExpZ;
    
    T density = 0;
    for(int k = 0; k < mix; ++k){
      density += pi(k) * dets(k) * pow(6.283185, -3) *  exp(-0.5 * kernel(k));
    }
    
    return log(density);
  }
};

double pLogDensMix(mat data, vec theta, mat mean, cube SigInv, vec dets, vec weights, int lagsA, int lagsD){
  int N = data.n_rows;
  int mix = weights.n_rows;
  // Constrained Positive
  double sigSqV = std::exp(theta(0)), sigSqD = std::exp(theta(1));
  // Evaluate log(p(theta)), start by evaluative the quadratic in the MVN exponents
  double prior = 0;
  for(int k = 0; k < mix; ++k){
    prior += weights(k) * pow(6.283185, -3) * dets(k) * exp(-0.5 * as_scalar((theta - mean.col(k)).t() * SigInv.slice(k) * (theta - mean.col(k))));
  }
  
  // Evaluate log likelihood
  double logLik = 0;
  int lags = std::max(lagsA, lagsD);
  for(int t = lags; t < N; ++t){
    double kernelA = data(t, 0);
    double kernelD = data(t, 1);
    for(int i = 1; i <= lagsA; ++i){
      kernelA -= data(t-i, 0) * theta(1 + i);
    }
    for(int i = 1; i <= lagsD; ++i){
      kernelD -= data(t-i, 1) * theta(1 + lagsA + i);
    }
    logLik += - 0.5 * (theta(0) + theta(1)) -
      std::pow(kernelA, 2) / (2 * sigSqV) - 
      std::pow(kernelD, 2) / (2 * sigSqD);
  }
  return std::log(prior) + logLik;
}

double pLogDensSingle(mat data, vec theta, vec mean, mat Linv, int lagsA, int lagsD){
  int N = data.n_rows;
  // Constrained Positive
  double sigSqV = std::exp(theta(0)), sigSqD = std::exp(theta(1));
  // Evaluate log(p(theta)), start by evaluative the quadratic in the MVN exponents
  
  // Evaluate log(p(theta))
  double prior = 0;
  mat sigInv = Linv.t() * Linv;
  prior = - 0.5 * as_scalar((theta - mean).t() * sigInv * (theta - mean));
  
  // Evaluate log likelihood
  // Evaluate log likelihood
  double logLik = 0;
  int lags = std::max(lagsA, lagsD);
  for(int t = lags; t < N; ++t){
    double kernelA = data(t, 0);
    double kernelD = data(t, 1);
    for(int i = 1; i <= lagsA; ++i){
      kernelA -= data(t-i, 0) * theta(1 + i);
    }
    for(int i = 1; i <= lagsD; ++i){
      kernelD -= data(t-i, 1) * theta(1 + lagsA + i);
    }
    logLik += - 0.5 * (theta(0) + theta(1)) -
      std::pow(kernelA, 2) / (2 * sigSqV) - 
      std::pow(kernelD, 2) / (2 * sigSqD);
  }
  return prior + logLik;
}

// These models are parameterised by the mean and log standard deviations
// [[Rcpp::export]]
Rcpp::List singlePriorMixApprox(mat data, Rcpp::NumericMatrix lambdaIn, vec theta, vec mean, mat Linv, int lagsA = 1, int lagsD = 1, int mix = 6){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1> grad(dim);
  
  mixQlogdens logQ(theta, mix);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(logQ, lambda, qEval, grad);
  
  double logp = pLogDensSingle(data, theta, mean, Linv, lagsA, lagsD);
  double elbo = logp - qEval;
  
  for(int i = 0; i < dim; ++i){
    grad(i) *= elbo;
  }
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = elbo);
}

// [[Rcpp::export]]
Rcpp::List mixPriorMixApprox(mat data, Rcpp::NumericMatrix lambdaIn, vec theta, mat mean, cube SigInv, vec dets, vec weights, int lagsA = 1, int lagsD = 1, int mix = 6){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1> grad(dim);
  
  mixQlogdens logQ(theta, mix);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(logQ, lambda, qEval, grad);
  
  double logp = pLogDensMix(data, theta, mean, SigInv, dets, weights, lagsA, lagsD);
  double elbo = logp - qEval;
  
  for(int i = 0; i < dim; ++i){
    grad(i) *= elbo;
  }
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = elbo);
}

// MCMC MH density calculation

// [[Rcpp::export]]
double nlogDensity (mat data, vec theta, vec mu, mat varInv, int lagsA = 1, int lagsD = 1){
  using namespace std;
  int T = data.n_rows;
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  double rho = 2 / (1 + exp(-theta(2))) - 1;
  
  int lags = max(lagsA, lagsD);
  for(int t = lags; t < T; ++t){
    
    double meanA = 0, meanD = 0;
    for(int i = 1; i <= lagsA; ++i){
      meanA += theta(2 + i) * data(t-i, 0);
    }
    for(int i = 1; i < lagsD; ++i){
      meanD += theta(2 + lagsA + i) * data(t - i, 1);
    }
    double z = pow(data(t, 0) - meanA, 2) / exp(theta(0))  +  pow(data(t, 1) - meanD, 2) / exp(theta(1))  -
      2 * rho * (data(t, 0) - meanA) * (data(t, 1) - meanD) / (sqrt(exp(theta(0))) * sqrt(exp(theta(1))));
    
    dens += - 0.5 * (theta(0) + theta(1) + log(1 - rho * rho))  -  z / (2 - 2 * rho * rho);
  }
  return dens;
}

// [[Rcpp::export]]
double nMixLogDens (mat data, vec theta, vec mu, mat varInv, vec weights, int lagsA = 1, int lagsD = 1){
  using namespace std;
  int T = data.n_rows;
  int dim = 2 + lagsA + lagsD;
  int mix = weights.n_elem;
  double dens = 0;
  for(int k = 0; k < mix; ++k){
    vec submean = mu.subvec(k*dim, (k+1)*dim-1);
    mat subinv = varInv.rows(k*dim, (k+1)*dim-1);
    dens += weights(k) * sqrt(det(subinv)) * exp(-0.5 * as_scalar((theta - submean).t() * subinv * (theta - submean)));
  }
  dens = std::log(dens);
  int lags = max(lagsA, lagsD);
  for(int t = lags; t < T; ++t){
    double meanA = 0, meanD = 0;
    for(int i = 1; i <= lagsA; ++i){
      meanA += theta(2 + i) * data(t-i, 0);
    }
    for(int i = 1; i < lagsD; ++i){
      meanD += theta(2 + lagsA + i) * data(t - i, 1);
    }
    double z = pow(data(t, 0) - meanA, 2) / exp(theta(0))  +  pow(data(t, 1) - meanD, 2) / exp(theta(1));
    
    dens += - 0.5 * (theta(0) + theta(1))  -  z / 2;
  }
  return dens;
}

// Forecast density evaluation

// [[Rcpp::export]]
cube evalVBDensIndep (mat data, vec mean, mat L, vec weights, int N, int H, mat grid, bool isMix, int dim, int K, int lagsA, int lagsD){
  cube densities (N, 2, H, fill::zeros);
  boost::math::normal_distribution<> Zdist(0, 1);
  
  for(int m = 0; m < 1000; ++m){
    // Set up lags of predicted means, starting with the true data
    int maxLags = std::max(lagsA, lagsD);
    vec alags = data(span(std::min(0, lagsD - lagsA), maxLags-1), 0), dlags = data(span(std::min(0, lagsA - lagsD), maxLags-1), 1);
    vec draws(dim);
    if(!isMix){
      // If not a mixture, draws = mu + L * eps
      draws = mean + L * randn<vec>(dim);
    } else {
      // If it is a mixture distribution choose one via inverse CDF weights sampling, then do mu_component + L_component * eps
      double u = randu<double>();
      int comp;
      for(int i = 0; i < K; ++i){
        if(u < weights(i)){
          comp = i;
          break;
        }
      }
      vec submean = mean.rows(comp*dim, (comp+1)*dim - 1);
      mat subL = L.rows(comp*dim, (comp+1)*dim-1);
      draws = submean + subL * randn<vec>(dim);
    }
    double sigmaSqA = std::exp(draws(0));
    double sigmaSqD = std::exp(draws(1));
    vec varA(lagsA, fill::zeros), varD(lagsD, fill::zeros);
    double afc, dfc;
    for(int h = 0; h < H; ++h){
      afc = dfc = 0;
      double currentVarA = sigmaSqA, currentVarD = sigmaSqD;
      for(int i = 1; i <= lagsA; ++i){
        afc += draws(1 + i) * alags(lagsA - i);
        currentVarA += pow(draws(1 + i), 2) * varA(lagsA - i);
      }
      for(int i = 1; i <= lagsD; ++i){
        dfc += draws(1 + lagsA + i) * dlags(lagsD - i);
        currentVarD += pow(draws(1 + lagsA + i), 2) * varD(lagsD - i);
      }
      double sdA = std::sqrt(currentVarA), sdD = std::sqrt(currentVarD);
      // Iterate densities over the grid
      for(int i = 0; i < N; ++i){
        densities(i, 0, h) += pdf(Zdist, (grid(i, 0) - afc) / sdA) / (1000 * sdA);
        densities(i, 1, h) += pdf(Zdist, (grid(i, 1) - dfc) / sdD) / (1000 * sdD);
      }
      // lag acceleration and angle
      for(int i = 0; i < lagsA - 1; ++i){
        alags(i) = alags(i+1);
        varA(i) = varA(i+1);
      }
      for(int i = 0; i < lagsD - 1; ++i){
        dlags(i) = dlags(i+1);
        varD(i) = varD(i+1);
      }
      alags(lagsA - 1) = afc;
      dlags(lagsD - 1) = dfc;
      varA(lagsA - 1) = currentVarA;
      varD(lagsD - 1) = currentVarD;
    }
  }
  return densities;    
}

// [[Rcpp::export]]
cube evalMCMCDensIndep (mat data, int N, int H, mat grid, mat MCMCdraws, int lagsA, int lagsD){
  cube densities (N, N, H, fill::zeros);
  boost::math::normal_distribution<> Zdist(0, 1);
  
  for(int m = 0; m < 1000; ++m){
    int maxLags = std::max(lagsA, lagsD);
    vec alags = data(span(std::min(0, lagsD - lagsA), maxLags-1), 0), dlags = data(span(std::min(0, lagsA - lagsD), maxLags-1), 1);
    vec draws = MCMCdraws.row(1000 + 4 * m).t();
    double sigmaSqA = std::exp(draws(0));
    double sigmaSqD = std::exp(draws(1));
    vec varA(lagsA, fill::zeros), varD(lagsD, fill::zeros);
    double afc, dfc;
    for(int h = 0; h < H; ++h){
      afc = dfc = 0;
      double currentVarA = sigmaSqA, currentVarD = sigmaSqD;
      for(int i = 1; i <= lagsA; ++i){
        afc += draws(1 + i) * alags(lagsA - i);
        currentVarA += pow(draws(1 + i), 2) * varA(lagsA - i);
      }
      for(int i = 1; i <= lagsD; ++i){
        dfc += draws(1 + lagsA + i) * dlags(lagsD - i);
        currentVarD += pow(draws(1 + lagsA + i), 2) * varD(lagsD - i);
      }
      double sdA = std::sqrt(currentVarA), sdD = std::sqrt(currentVarD);
      // Iterate densities over the grid
      for(int i = 0; i < N; ++i){
        densities(i, 0, h) += pdf(Zdist, (grid(i, 0) - afc) / sdA) / (1000 * sdA);
        densities(i, 1, h) += pdf(Zdist, (grid(i, 1) - dfc) / sdD) / (1000 * sdD);
      }
      // lag acceleration and angle
      for(int i = 0; i < lagsA - 1; ++i){
        alags(i) = alags(i+1);
        varA(i) = varA(i+1);
      }
      for(int i = 0; i < lagsD - 1; ++i){
        dlags(i) = dlags(i+1);
        varD(i) = varD(i+1);
      }
      alags(lagsA - 1) = afc;
      dlags(lagsD - 1) = dfc;
      varA(lagsA - 1) = currentVarA;
      varD(lagsD - 1) = currentVarD;
    }
  }
  return densities;
}

