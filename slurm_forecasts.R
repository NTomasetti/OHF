rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(repenv)
set.seed(1000 + i)

library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
library(RcppEigen, lib.loc = 'packages')
library(rstan, lib.loc = 'packages')
library(ltsa, lib.loc = 'packages')
source('slurm_RFuns.R')
sourceCpp('slurm_cppFuns.cpp')

ids <- readRDS('idTest.RDS')
datatest <- readRDS('dataTest.RDS')

data <- datatest[[i]]
id <- ids[i]

homogDraws <- readRDS('homogMCMCN2000.RDS')$draws[seq(2501, 7500, 5),]
CHprior <- readRDS('CHFit.RDS')

H <- 30
S <- 10
minT <- S
maxT <- 400
sSeq <- seq(minT, maxT, S)
results <- data.frame()
lagsA <- 4
lagsD <- 4
K <- mix <- 6
dim <- 2 + lagsA + lagsD

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

starting <- list(matrix(c(prior[[2]]$mean[1:10], diag(0.5, 10)), ncol = 1),
                 matrix(c(prior[[2]]$mean, rep(log(0.5), 10*K), rep(1, K)), ncol = 1))

# Set forcast supports
aLower <- min(data[,1])
aUpper <- max(data[,1])
dLower <- min(data[,2])
dUpper <- max(data[,2])
asup <- seq(aLower - 0.1, aUpper + 0.1, length.out = 200)
dsup <- seq(dLower - 0.1, dUpper + 0.1, length.out = 200)
grid <- cbind(asup, dsup)

# Incrementally add data to VB fits
for(s in seq_along(sSeq)){
  if((sSeq[s] + H) > nrow(data)){
    break
  }

  # Offline (Standard VB) model fits
  fitOffline <- fitVB(data[1:sSeq[s], 1:2], prior, starting, dim = 10, mix = K)

  if(s == 1){
    fitOnline <- fitOffline
  } else {
    fitOnline <- fitUVB(data[(sSeq[s-1]-(max(lagsA, lagsD) - 1)):sSeq[s], 1:2], fitOnline, starting, dim = 10, mix = K)
  }

  # Get MCMC posteriors
  MCMCDraw <- list()
  for(k in 1:2){
    MCMCDraw[[k]] <- singleMCMCallMH(data[1:sSeq[s],1:2], 15000, homogDraws[1, ], hyper[[k]],
                                 stepsize = 0.05, mix = (k == 2))$draws[seq(10001, 15000, 5),]
  }
  
  # Evaluate predictive densities
  SVB <- lapply(1:4, function(x) VBDens(data = data[(sSeq[s]-(max(lagsA, lagsD) - 1)):(sSeq[s] + H), ],
                                        fit = fitOffline[[x]],
                                        grid = grid, 
                                        H = H,
                                        dim = 10,
                                        K = K,
                                        isMix = (x %% 2 == 0),
                                        lagsA = lagsA,
                                        lagsD = lagsD))
 
  if(s == 1){
    UVB <- SVB
  } else {
    UVB <- lapply(1:4, function(x) VBDens(data = data[(sSeq[s]-(max(lagsA, lagsD) - 1)):(sSeq[s] + H), ],
                                          fit = fitOnline[[x]],
                                          grid = grid, 
                                          H = H,
                                          dim = 10,
                                          K = K,
                                          isMix = (x %% 2 == 0),
                                          lagsA = lagsA,
                                          lagsD = lagsD))
  }
  
  MCMC <- lapply(MCMCDraw, function(x) MCMCDens(data = data[(sSeq[s]-(max(lagsA, lagsD)-1)):(sSeq[s] + H), ],
                                            H = H,
                                            grid = grid,
                                            MCMCdraws = x,
                                            lagsA = lagsA,
                                            lagsD = lagsD))
  
  # Grab logscores etc. for heterogenous models
  results <- rbind(results,
                   data.frame(logscore = c(SVB[[1]]$logscore, SVB[[2]]$logscore, SVB[[3]]$logscore, SVB[[4]]$logscore,
                                           UVB[[1]]$logscore, UVB[[2]]$logscore, UVB[[3]]$logscore, UVB[[4]]$logscore,
                                           MCMC[[1]]$logscore, MCMC[[2]]$logscore),
                              dist = c(SVB[[1]]$mapDist, SVB[[2]]$mapDist, SVB[[3]]$mapDist, SVB[[4]]$mapDist,
                                       UVB[[1]]$mapDist, UVB[[2]]$mapDist, UVB[[3]]$mapDist, UVB[[4]]$mapDist,
                                       MCMC[[1]]$mapDist, MCMC[[2]]$mapDist),
                              h = rep(1:H, 10),
                              model = rep(c('IH', 'IH', 'CH', 'CH', 'IH', 'IH', 'CH', 'CH', 'IH', 'CH'), rep(H, 10)),
                              inference = rep(c('SVB-Single', 'SVB-Mixture', 'SVB-Single', 'SVB-Mixture', 'UVB-Single', 
                                                'UVB-Mixture', 'UVB-Single', 'UVB-Mixture', 'MCMC', 'MCMC'), rep(H, 10)),
                              earthmover = c(SVB[[1]]$em, SVB[[2]]$em, SVB[[3]]$em, SVB[[4]]$em, 
                                            UVB[[1]]$em, UVB[[2]]$em, UVB[[3]]$em, UVB[[4]]$em, rep(NA, 2*H)),
                              S = sSeq[s],
                              id = id)
                   )
  # Homogenous model results
  homogResults <-  MCMCDens(data = data[(sSeq[s]-(max(lagsA, lagsD) - 1)):(sSeq[s] + H), ],
                            H = H,
                            grid = grid,
                            MCMCdraws = homogDraws,
                            lagsA = lagsA,
                            lagsD = lagsD)
  
  results <- rbind(results,
                   data.frame(logscore = homogResults$logscore,
                              dist = homogResults$mapDist,
                              h = 1:H,
                              model = 'Homogenous',
                              inference = 'MCMC',
                              earthmover = NA,
                              S = sSeq[s],
                              id = id))
  
  # Naive Model Forecasts
  dConst <- data[sSeq[s], 2]
  aConst <- data[sSeq[s] +1, 3] - data[sSeq[s] , 3]
  aAvg <- 0.1 * (data[sSeq[s] + 1, 3] - data[sSeq[s] -9, 3])
  dAvg <- mean(data[-9:0 + sSeq[s], 2])
  # Model 1: Past A, Const D
  # Model 2: Past A, Zero D
  # Model 3: Past A, Mean D
  # Model 4-6: Mean A, same order of D
  # Model 7-9: Const V (Zero A), same order of D
  x0 <- rep(data[sSeq[s], 4], 9)
  y0 <- rep(data[sSeq[s], 5], 9)
  v0 <- rep(data[sSeq[s], 3], 9) 
  
  
  for(h in 1:H){
    v0 <- v0 + c(rep(aConst, 3), rep(aAvg, 3), rep(0, 3))
    x0 <- x0 + v0 * rep(c(cos(dConst + pi/2), cos(pi/2), cos(pi/2 + dAvg)), 3)
    y0 <- y0 + v0 * rep(c(sin(dConst + pi/2), sin(pi/2), sin(pi/2 + dAvg)), 3)
    dist <- sqrt((x0 - data[sSeq[s] + h, 4])^2 + (y0 - data[sSeq[s] + h, 5])^2)
    results <- rbind(results,
                     data.frame(logscore = NA,
                                dist = dist,
                                h = h,
                                model = paste('Naive', 1:9),
                                inference = NA,
                                earthmover = NA,
                                S = sSeq[s],
                                id = id))
  }
}

write.csv(results, paste0('eval/car', id, '.csv'), row.names=FALSE)

