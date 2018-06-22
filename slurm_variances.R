rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
j <- as.numeric(repenv)
set.seed(1000 + j)

library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
library(RcppEigen, lib.loc = 'packages')
library(rstan, lib.loc = 'packages')
source('slurm_RFuns.R')
sourceCpp('slurm_cppFuns.cpp')

ids <- readRDS('idTest.RDS')
datatest <- readRDS('dataTest.RDS')

results <- data.frame()

homogDraws <- readRDS('homogMCMCN2000.RDS')$draws[floor(seq(5001, 7500, length.out = 500)),]
CHprior <- readRDS('CHFit.RDS')

T <- 400
lagsA <- 4
lagsD <- 4
K <- mix <- 6
dim <- 2 + lagsA + lagsD
model <- c('IH', 'CH')
# VB Prior Distributions: 1) IH, 2) CH
prior <- list()
prior[[1]] <- c(-5, -5, rep(0, 8), c(chol(diag(10, 10))))
prior[[2]] <- CHprior

# MCMC Hyper parameters
hyper <- list()
hyper[[1]] <- list()
hyper[[1]]$mean <- c(-5, -5, rep(0, 8))
hyper[[1]]$varInv <- solve(diag(10, 10))
hyper[[2]] <- list(mean = prior[[2]]$mean, varInv = prior[[2]]$varInv, weights = prior[[2]]$weights)

for(i in 1:100){
  
  data <- datatest[[(j-1)*100 + i]]
  id <- ids[(j-1)*100 + i]
  
  # Get MCMC posteriors
  MCMCDraw <- list()
  for(k in 1:2){
    MCMCDraw[[k]] <- singleMCMCallMH(data = data[1:min(nrow(data), T), 1:2],
                                     reps = 15000,
                                     draw = homogDraws[1, ],
                                     hyper = hyper[[k]],
                                     stepsize = 0.05,
                                     mix = (k == 2))$draws[floor(seq(10001, 15000, length.out = 500)),]
   
    var <- colMeans(exp(MCMCDraw[[k]][,1:2]))
    results <- rbind(results,
                     data.frame(id = id,
                                model = model[k],
                                variance = var,
                                parameter = c('a', 'delta')))
     
  }
}



write.csv(results, paste0('variances/group', j, '.csv'), row.names=FALSE)

