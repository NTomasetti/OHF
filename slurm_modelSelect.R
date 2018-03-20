rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(repenv)
set.seed(1000 + i)

library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
library(tidyverse, lib.loc = 'packages')
source('MCMCFuns.R')

dataset <- readRDS('TrainData.RDS')
data <- datatest[[i]]

# Loop through VAR: TRUE / FALSE and up to six lags for both variables to
# 1) Fit independent MCMC for data from each vehicle for 1:(T-100)
# 2) Calculate posterior means for each element of theta
# 3) Conduct one step ahead forecasts from T-99:T
# 4) Calculate squared error
reps <- 20000
thin <- 10
var <- c(TRUE, FALSE)

# Save both the squared forecast error and posterior means.
error <- NULL
for(k in 1:2){
  for(lags in 1:6){
    T <- nrow(data)
    # Var = TRUE and FALSE have different parameters
    if(var[k]){
      vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'rho', paste0('phi[', 1:(2*lags), ']'), paste0('gamma[', 1:(2*lags), ']'))
      draw <- c(-5, -5, rep(0, 1 + 4 * lags))
      hyper <- list(mean = c(-5, -5, rep(0, 1 + 4 * lags)), varInv = solve(diag(c(10, 10, 2.5, rep(10, 4 * lags)))))
    } else {
      vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'rho', paste0('phi[', 1:lags, ']'), paste0('gamma[', 1:lags, ']'))
      draw <- c(-5, -5, rep(0, 1 + 2 * lags))
      hyper <- list(mean = c(-5, -5, rep(0, 1 + 2 * lags)), varInv = solve(diag(c(10, 10, 2.5, rep(10, 2 * lags)))))
    }
    MCMC <- as.data.frame(indepMCMC(data[1:(T-100),], reps, draw, hyper, lagsA = lags, lagsD = lags, thin = thin, VAR = var[k])$draws)
      
    # Transform the variance matrix from R^3 to sigma^2 and rho
    colnames(MCMC) <- vars
    MCMC %>%
      mutate(iter = seq_along(rho) * thin) %>%
      filter(iter > reps/2) %>%
      mutate(`sigma[epsilon]^{2}` = exp(`sigma[epsilon]^{2}`),
             `sigma[eta]^{2}` = exp(`sigma[eta]^{2}`),
             rho = 2 / (1 + exp(-rho)) - 1) %>%
      select(-iter) %>%
      colMeans() -> means

    for(t in (T-100):(T-10)){
      lagA <- data[0:(1 - lags) + t, 1]
      lagD <- data[0:(1 - lags) + t, 2]
    
      for(h in 1:10){
        a <- d <- 0
        for(j in 1:lags){
          if(var[k]){
            a <- a  +  lagA[j] * means[3 + j]  +  abs(lagD[j]) * sign(lagA[j]) * means[3 + lags + j]
            d <- d  +  lagD[j] * means[3 + 2 * lags + j]  +  abs(lagA[j]) * sign(lagD[j]) * means[3 + 3 * lags + j]
          } else {
            a <- a  +  lagA[j] * means[3 + j] 
            d <- d  +  lagD[j] * means[3 + lags + j]
          }
        }
        lagA <- c(a, lagA[1:(lags - 1)])
        lagD <- c(d, lagD[1:(lags - 1)])
      }
      
      error <- rbind(error, data.frame(a = (a - data[t+10, 1])^2,
                                       d = (d - data[t+10, 2])^2,
                                       lags = lags,
                                       VAR = var[k],
                                       ID = i))
    }
  }
}
# Forecasts assuming a constant mean of 0, pi/2 (note: d is already delta - pi/2) and no estimated lags.
T <- nrow(data)
for(t in (T-100):(T-10)){
  error <- rbind(error, data.frame(a = data[t+10, 1]^2,
                                   d = data[t+10, 2]^2,
                                   lags = 0,
                                   VAR = c(TRUE, FALSE),
                                   ID = i))
}

write.csv(error, paste0('modelSelection/car', i, '.csv'), row.names = FALSE)