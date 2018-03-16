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

//TODO: Change likelihood to dependent structure (and possibly another model entirely)

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
  const int lags;
  SPSA(const mat& dataIn, const vec& epsIn, const vec& meanIn, const mat& LinvIn, const int& lagsIn) :
    data(dataIn), epsilon(epsIn), mean(meanIn), Linv(LinvIn), lags(lagsIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    int N = data.n_rows;
    int dim = 2 + 2 * lags;
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
    for(int t = lags; t < N; ++t){
      loglikKernel(t, 0) = data(t, 0);
      loglikKernel(t, 1) = data(t, 1);
      for(int i = 1; i <= lags; ++i){
        loglikKernel(t, 0) -= data(t-i, 0) * theta(1 + i);
        loglikKernel(t, 1) -= data(t-i, 1) * theta(1 + lags + i);
      }
      logLik += - 0.5 * theta(0) - 0.5 * theta(1) - 
        pow(loglikKernel(t, 0), 2) / (2 * sigSqV) -
        pow(loglikKernel(t, 1), 2) / (2 * sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List singlePriorSingleApprox(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec mean, mat Linv, int lags){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1>  grad(6 + 10 * lags + 4 * lags * lags);
  // Autodiff
  SPSA p(data, epsilon, mean, Linv, lags);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

// Mixture Prior / Single Approx

struct MPSA {
  const mat data;
  const vec epsilon;
  const vec mean;
  const cube sigInv;
  const vec dets;
  const vec weights;
  const int lags;
  MPSA(const mat& dataIn, const vec& epsIn, const vec& meanIn, const cube& sigInvIn, const vec& detsIn, const vec& weightsIn, const int& lagsIn) :
    data(dataIn), epsilon(epsIn), mean(meanIn), sigInv(sigInvIn), dets(detsIn), weights(weightsIn), lags(lagsIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    int N = data.n_rows;
    int dim = 2 + 2 * lags;
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
    for(int k = 0; k < 6; ++k){
      prior += weights(k) * pow(6.283185, -3) * dets(k) * exp(-0.5 * as_scalar((theta - mean.col(k)).t() * sigInv.slice(k) * (theta - mean.col(k))));
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
    for(int t = lags; t < N; ++t){
      loglikKernel(t, 0) = data(t, 0);
      loglikKernel(t, 1) = data(t, 1);
      for(int i = 1; i <= lags; ++i){
        loglikKernel(t, 0) -= data(t-i, 0) * theta(1 + i);
        loglikKernel(t, 1) -= data(t-i, 1) * theta(1 + lags + i);
      }
      logLik += - 0.5 * theta(0) - 0.5 * theta(1) - 
        pow(loglikKernel(t, 0), 2) / (2 * sigSqV) -
        pow(loglikKernel(t, 1), 2) / (2 * sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List mixPriorSingleApprox(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec mean, cube sigInv, vec dets, vec weights, int lags){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1>  grad(6 + 10 * lags + 4 * lags * lags);
  // Autodiff
  MPSA p(data, epsilon, mean, sigInv, dets, weights, lags);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

// Pieces for the mixture approximations are shared between both prior specifications

struct mixQlogdens {
  const vec theta;
  mixQlogdens(const vec& thetaIn) :
    theta(thetaIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt;
    
    Matrix<T, Dynamic, 1> dets(6);
    for(int k = 0; k < 6; ++k){
      dets(k) = exp(lambda(6*k + 36));
      for(int i = 1; i < 6; ++i){
        dets(k) *= exp(lambda(36 + 6*k + i));
      }
      dets(k) = 1.0 / dets(k);
    }
    
    Matrix<T, Dynamic, 1> kernel(6);
    kernel.fill(0);
    
    for(int k = 0; k < 6; ++k){
      for(int i = 0; i < 6; ++i){
        kernel(k) += pow((theta(i) - lambda(k*6 + i)) / exp(lambda(36 + 6*k + i)), 2);
      }
    }
    
    Matrix<T, Dynamic, 1> pi(6);
    T sumExpZ = 0;
    for(int k = 0; k < 6; ++k){
      pi(k) = exp(lambda(72 + k));
      sumExpZ += pi(k);
    }
    pi /= sumExpZ;
    
    T density = 0;
    for(int k = 0; k < 6; ++k){
      density += pi(k) * dets(k) * pow(6.283185, -3) *  exp(-0.5 * kernel(k));
    }
    
    return log(density);
  }
};

double pLogDensMix(mat data, vec theta, mat mean, cube SigInv, vec dets, vec weights){
  int N = data.n_rows;
  // Constrained Positive
  double sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
  // Evaluate log(p(theta)), start by evaluative the quadratic in the MVN exponents
  double prior = 0;
  for(int k = 0; k < 6; ++k){
    prior += weights(k) * pow(6.283185, -3) * dets(k) * exp(-0.5 * as_scalar((theta - mean.col(k)).t() * SigInv.slice(k) * (theta - mean.col(k))));
  }
  
  // Evaluate log likelihood
  double logLik = 0;
  for(int t = 2; t < N; ++t){
    logLik += - 0.5 * (theta(0) + theta(1)) - pow(data(t, 0) - theta(2) * data(t-1, 0) - theta(3) * data(t-2, 0), 2) / (2 * sigSqV) -
      pow(data(t, 1) - theta(4) * data(t-1, 1) - theta(5) * data(t-2, 1), 2) / (2 * sigSqD);
  }
  return log(prior) + logLik;
}

double pLogDensSingle(mat data, vec theta, vec mean, mat Linv){
  int N = data.n_rows;
  // Constrained Positive
  double sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
  // Evaluate log(p(theta)), start by evaluative the quadratic in the MVN exponents
  
  // Evaluate log(p(theta))
  double prior = 0;
  mat sigInv = Linv.t() * Linv;
  prior = - 0.5 * as_scalar((theta - mean).t() * sigInv * (theta - mean));
  
  // Evaluate log likelihood
  double logLik = 0;
  for(int t = 2; t < N; ++t){
    logLik += - 0.5 * (theta(0) + theta(1)) - pow(data(t, 0) - theta(2) * data(t-1, 0) - theta(3) * data(t-2, 0), 2) / (2 * sigSqV) -
      pow(data(t, 1) - theta(4) * data(t-1, 1) - theta(5) * data(t-2, 1), 2) / (2 * sigSqD);
  }
  return prior + logLik;
}

// These models are parameterised by the mean and log standard deviations
// [[Rcpp::export]]
Rcpp::List mixPriorMixApprox(mat data, Rcpp::NumericMatrix lambdaIn, vec theta, mat mean, cube SigInv, vec dets, vec weights){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  Matrix<double, Dynamic, 1> grad(13*6);
  
  mixQlogdens logQ(theta);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(logQ, lambda, qEval, grad);
  
  double logp = pLogDensMix(data, theta, mean, SigInv, dets, weights);
  double elbo = logp - qEval;
  
  for(int i = 0; i < 13*6; ++i){
    grad(i) *= elbo;
  }
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = elbo);
}

// [[Rcpp::export]]
Rcpp::List singlePriorMixApprox(mat data, Rcpp::NumericMatrix lambdaIn, vec theta, vec mean, mat Linv){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  Matrix<double, Dynamic, 1> grad(13*6);
  
  mixQlogdens logQ(theta);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(logQ, lambda, qEval, grad);
  
  double logp = pLogDensSingle(data, theta, mean, SigInv, dets, weights);
  double elbo = logp - qEval;
  
  for(int i = 0; i < 13*6; ++i){
    grad(i) *= elbo;
  }
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = elbo);
}

// MCMC MH density calculation

// [[Rcpp::export]]
double nlogDensity (mat data, vec theta, vec mu, mat varInv, int lags = 1){
  int T = data.n_rows;
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  double rho = 2 / (1 + exp(-theta(2))) - 1;
  
  for(int t = lags; t < T; ++t){
    
    double meanA = 0, meanD = 0;
    for(int i = 1; i <= lags; ++i){
      meanA += theta(2 + i) * data(t-i, 0);
      meanD += theta(2 + lags + i) * data(t - i, 1);
    }
    double z = pow(data(t, 0) - meanA, 2) / exp(theta(0))  +  pow(data(t, 1) - meanD, 2) / exp(theta(1))  -
      2 * rho * (data(t, 0) - meanA) * (data(t, 1) - meanD) / (sqrt(exp(theta(0))) * sqrt(exp(theta(1))));
    
    dens += - 0.5 * (theta(0) + theta(1) + log(1 - rho * rho))  -  z / (2 - 2 * rho * rho);
  }
  return dens;
}

// [[Rcpp::export]]
double nMixLogDens (mat data, vec theta, vec mu, mat varInv, vec weights, int lags = 1){
  int T = data.n_rows;
  int K = weights.n_elem;
  double dens = 0;
  for(int k = 0; k < K; ++k){
    vec submean = mu.subvec(k*6, (k+1)*6-1);
    mat subinv = varInv.rows(k*6, (k+1)*6-1);
    dens += weights(k) * std::sqrt(det(subinv)) * std::exp(-0.5 * as_scalar((theta - submean).t() * subinv * (theta - submean)));
  }
  dens = std::log(dens);
  for(int t = lags; t < T; ++t){
    
    double meanA = 0, meanD = 0;
    for(int i = 1; i <= lags; ++i){
      meanA += theta(2 + i) * data(t-i, 0);
      meanD += theta(2 + lags + i) * data(t - i, 1);
    }
    double z = pow(data(t, 0) - meanA, 2) / exp(theta(0))  +  pow(data(t, 1) - meanD, 2) / exp(theta(1))  -
      2 * rho * (data(t, 0) - meanA) * (data(t, 1) - meanD) / (sqrt(exp(theta(0))) * sqrt(exp(theta(1))));
    
    dens += - 0.5 * (theta(0) + theta(1) + log(1 - rho * rho))  -  z / (2 - 2 * rho * rho);
  }
  return dens;
}

// Forecast density evaluation

// [[Rcpp::export]]
cube evalMCMCDens (mat data, int N, int H, mat grid, mat MCMCdraws){
  cube densities (N, N, H, fill::zeros);
  boost::math::normal_distribution<> Zdist(0, 1);
  
  for(int m = 0; m < 1000; ++m){
    double afc1 = data(1, 0), afc2 = data(0, 0), dfc1 = data(1, 1), dfc2 = data(0, 1), afc, dfc, vfc;
    vec draws = MCMCdraws.row(1000 + 4 * m).t();
    double sigV = std::sqrt(std::exp(draws(0)));
    double sigD = std::sqrt(std::exp(draws(1)));
    for(int h = 0; h < H; ++h){
      afc = afc1 * draws(2) + afc2 * draws(3);
      dfc = dfc1 * draws(4) + dfc2 * draws(5);
      vfc = afc + data(2 + h, 2);
      // Iterate densities over the grid
      for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
          // Transform delta x / delta y to velocity and angle
          double v = std::sqrt(
            std::pow(grid(i*N + j, 0), 2) +
              std::pow(grid(i*N + j, 1), 2)
          );
          double del = std::atan2(
            grid(i*N + j, 1),
            grid(i*N + j, 0)
          ) - 1.570796;
          // p(x, y) = p(v) * p(d) / sqrt(dx^2 + dy^2)
          densities(i, j, h) += (pdf(Zdist, (v - vfc) / sigV) / sigV) *
            (pdf(Zdist, (del - dfc) / sigD) / sigD) / (v * 1000);
        }
      }
      afc2 = afc1;
      afc1 = afc;
      dfc2 = dfc1;
      dfc1 = dfc;
    }
  }
  return densities;
}

// [[Rcpp::export]]
cube evalVBDens (mat data, vec mean, mat L, vec weights, int N, int H, mat grid, bool mix){
  cube densities (N, N, H, fill::zeros);
  boost::math::normal_distribution<> Zdist(0, 1);
  
  for(int m = 0; m < 1000; ++m){
    // Set up lags of predicted means, starting with the true data
    double afc1 = data(1, 0), afc2 = data(0, 0), dfc1 = data(1, 1), dfc2 = data(0, 1), afc, dfc, vfc;
    vec draws(6);
    if(!mix){
      // If not a mixture, draws = mu + L * eps
      draws = mean + L * randn<vec>(6);
    } else {
      // If it is a mixture distribution choose one via inverse CDF weights sampling, then do mu_component + L_component * eps
      double u = randu<double>();
      int comp;
      for(int i = 0; i < 6; ++i){
        if(u < weights(i)){
          comp = i;
          break;
        }
      }
      vec submean = mean.rows(comp*6, (comp+1)*6 - 1);
      mat subL = L.rows(comp*6, (comp+1)*6-1);
      draws = submean + subL * randn<vec>(6);
    }
    double sigV = std::sqrt(std::exp(draws(0)));
    double sigD = std::sqrt(std::exp(draws(1)));
    for(int h = 0; h < H; ++h){
      // Step afc and dfc means forward, vfc = afc + lagged velocity
      afc = afc1 * draws(2) + afc2 * draws(3);
      dfc = dfc1 * draws(4) + dfc2 * draws(5);
      vfc = afc + data(2 + h, 2);
      // Iterate densities over the grid
      for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
          // Transform delta x / delta y to velocity and angle
          double v = std::sqrt(
            std::pow(grid(i*N + j, 0), 2) +
              std::pow(grid(i*N + j, 1), 2)
          );
          double del = std::atan2(
            grid(i*N + j, 1),
            grid(i*N + j, 0)
          ) - 1.570796;
          // p(x, y) = p(v) * p(d) / sqrt(dx^2 + dy^2)
          densities(i, j, h) += (pdf(Zdist, (v - vfc) / sigV) / sigV) *
            (pdf(Zdist, (del - dfc) / sigD) / sigD) / (v * 1000);
        }
      }
      // lag acceleration and angle
      afc2 = afc1;
      afc1 = afc;
      dfc2 = dfc1;
      dfc1 = dfc;
    }
  }
  
  return densities;    
}

// [[Rcpp::export]]
cube evalVBDensIndep (mat data, vec mean, mat L, vec weights, int N, int H, mat grid, bool mix){
  cube densities (N, 2, H, fill::zeros);
  boost::math::normal_distribution<> Zdist(0, 1);
  
  for(int m = 0; m < 1000; ++m){
    // Set up lags of predicted means, starting with the true data
    double afc1 = data(1, 0), afc2 = data(0, 0), dfc1 = data(1, 1), dfc2 = data(0, 1), afc, dfc, vfc;
    vec draws(6);
    if(!mix){
      // If not a mixture, draws = mu + L * eps
      draws = mean + L * randn<vec>(6);
    } else {
      // If it is a mixture distribution choose one via inverse CDF weights sampling, then do mu_component + L_component * eps
      double u = randu<double>();
      int comp;
      for(int i = 0; i < 6; ++i){
        if(u < weights(i)){
          comp = i;
          break;
        }
      }
      vec submean = mean.rows(comp*6, (comp+1)*6 - 1);
      mat subL = L.rows(comp*6, (comp+1)*6-1);
      draws = submean + subL * randn<vec>(6);
    }
    double sigV = std::sqrt(std::exp(draws(0)));
    double sigD = std::sqrt(std::exp(draws(1)));
    for(int h = 0; h < H; ++h){
      // Step afc and dfc means forward, vfc = afc + lagged velocity
      afc = afc1 * draws(2) + afc2 * draws(3);
      dfc = dfc1 * draws(4) + dfc2 * draws(5);
      // Iterate densities
      for(int i = 0; i < N; ++i){
        densities(i, 0, h) += pdf(Zdist, (grid(i, 0) - afc) / sigV) / (1000 * sigV);
        densities(i, 1, h) += pdf(Zdist, (grid(i, 1) - dfc) / sigD) / (1000 * sigD);
      }
      // lag acceleration and angle
      afc2 = afc1;
      afc1 = afc;
      dfc2 = dfc1;
      dfc1 = dfc;
    }
  }
  return densities;    
}

// [[Rcpp::export]]
cube evalMCMCDensIndep (mat data, int N, int H, mat grid, mat MCMCdraws){
  cube densities (N, N, H, fill::zeros);
  boost::math::normal_distribution<> Zdist(0, 1);
  
  for(int m = 0; m < 1000; ++m){
    double afc1 = data(1, 0), afc2 = data(0, 0), dfc1 = data(1, 1), dfc2 = data(0, 1), afc, dfc, vfc;
    vec draws = MCMCdraws.row(1000 + 4 * m).t();
    double sigV = std::sqrt(std::exp(draws(0)));
    double sigD = std::sqrt(std::exp(draws(1)));
    for(int h = 0; h < H; ++h){
      afc = afc1 * draws(2) + afc2 * draws(3);
      dfc = dfc1 * draws(4) + dfc2 * draws(5);
      // Iterate densities over the grid
      for(int i = 0; i < N; ++i){
        densities(i, 0, h) += pdf(Zdist, (grid(i, 0) - afc) / sigV) / (1000 * sigV);
        densities(i, 1, h) += pdf(Zdist, (grid(i, 1) - dfc) / sigD) / (1000 * sigD);
      }
      afc2 = afc1;
      afc1 = afc;
      dfc2 = dfc1;
      dfc1 = dfc;
    }
  }
  return densities;
}

