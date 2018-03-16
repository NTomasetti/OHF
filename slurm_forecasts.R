rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(repenv)
set.seed(1000 + i)

library(Rcpp)#, lib.loc = 'packages')
library(RcppArmadillo)#, lib.loc = 'packages')
library(RcppEigen)#, lib.loc = 'packages')
library(rstan)#, lib.loc = 'packages')
source('slurmRFuns.R')
sourceCpp('slurmCppFuns.cpp')

ids <- readRDS('testID.RDS')

datatest <- readRDS('TestData.RDS')
data <- datatest[[i]]
id <- ids[i]

homogDraws <- readRDS('noHierN2000.RDS')$draws
prior <- readRDS('prior.RDS')

H <- 30
S <- 10
maxT <- 450
sSeq <- seq(100, maxT, S)
laneSize <- 3.5
results <- data.frame()


# MCMC Hyper parameters
hyper <- list()
hyper[[1]] <- list()
hyper[[1]]$mean <- prior[[1]][1:6]
uinv <- solve(matrix(prior[[1]][7:42], 6))
hyper[[1]]$varInv <- t(uinv) %*% uinv

# Mixture Prior - The prior.RDS object has mean / inverse u however this is converted to mean / log(sd) for the diagonal mixture model in VB
mean <- prior[[2]][1:36]
priorDiag <- mean
varinv <- NULL
for(k in 1:6){
  uinv <- matrix(prior[[2]][k*36 + 1:36], 6)
  vari <- t(uinv) %*% uinv
  logsd <- log(diag(solve(vari)))
  varinv <- rbind(varinv, vari)
  priorDiag <- c(priorDiag, logsd)
}
weights <- prior[[2]][6*7*6 + 1:6]
priorDiag <- c(priorDiag, weights)
weights <- exp(weights) / sum(exp(weights))
hyper[[2]] <- list(mean = mean, varInv = varinv, weights = weights)
prior[[2]] <- matrix(priorDiag, ncol = 1)

starting <- list(matrix(c(-5, -5, rep(0, 4), c(chol(diag(0.5, 6)))), ncol = 1),
                 prior[[2]])

# Set forcast supports
aLower <- min(data[,1])
if(aLower < 0){
  aLower <- 2 * aLower
} else {
  aLower <- aLower - 0.5
}
dLower <- min(data[,2])
if(dLower < 0){
  dLower <- 2 * dLower
} else {
  dLower <- 0.5 * dLower
}

asup <- seq(aLower, 2*max(data[,1]), length.out = 300)
dsup <- seq(dLower, 2*max(data[,2]), length.out= 300)
grid <- cbind(asup, dsup)

# Incrementally add data to VB fits
for(s in seq_along(sSeq)){
  if((sSeq[s] + H) > nrow(data)){
    break
  }

  # Offline (Standard VB) model fits
  fitOffline <- fitCarMods(data[1:sSeq[s], 1:2], prior, starting)

  if(s == 1){
    fitOnline <- fitOffline
  } else {
    fitOnline <- fitCarMods(data[(sSeq[s-1]-1):sSeq[s], 1:2], prior, starting)
  }

  # Get MCMC posteriors
  MCMC <- list()
  for(k in 1:2){
    MCMC[[k]] <- singleMCMCallMH(data[1:sSeq[s],1:2], 5000, c(-5, -5, 0, 0, 0, 0), hyper[[k]],
                                 stepsize = 0.05, mix = (k == 2))$draws
  }
  
  # Evaluate predictive densities
  SVB <- lapply(1:4, function(x) VBDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                        fit = fitOffline[[x]],
                                        grid = grid, 
                                        H = H,
                                        mix = (x %% 2 == 0),
                                        laneSize = laneSize))
 
  if(s == 1){
    UVB <- resultsOffline
  } else {
    UVB <- lapply(1:4, function(x) VBDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                          fit = fitOnline[[x]],
                                          grid = grid, 
                                          H = H,
                                          mix = (x %% 2 == 0),
                                          laneSize = laneSize))
    
   
  }
  MCMC <- lapply(MCMC, function(x) MCMCDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                            N = 300,
                                            H = H,
                                            grid = grid,
                                            MCMCdraws = x,
                                            laneSize = laneSize))
  
  # Grab logscores etc. for heterogenous models
  results <- rbind(results,
                   data.frame(logscore = c(SVB[[1]]$logscore, SVB[[2]]$logscore, SVB[[3]]$logscore, SVB[[4]]$logscore,
                                           UVB[[1]]$logscore, UVB[[2]]$logscore, UVB[[3]]$logscore, UVB[[4]]$logscore,
                                           MCMC[[1]]$logscore, MCMC[[2]]$logscore),
                              xCDF = c(SVB[[1]]$xCDF, SVB[[2]]$xCDF, SVB[[3]]$xCDF, SVB[[4]]$xCDF, 
                                       UVB[[1]]$xCDF, UVB[[2]]$xCDF, UVB[[3]]$xCDF, UVB[[4]]$xCDF,
                                       MCMC[[1]]$xCDF, MCMC[[2]]$xCDF),
                              dist = c(SVB[[1]]$mapDist, SVB[[2]]$mapDist, SVB[[3]]$mapDist, SVB[[4]]$mapDist,
                                       UVB[[1]]$mapDist, UVB[[2]]$mapDist, UVB[[3]]$mapDist, UVB[[4]]$mapDist,
                                       MCMC[[1]]$mapDist, MCMC[[2]]$mapDist),
                              probLeft = c(SVB[[1]]$probLeft, SVB[[2]]$probLeft, SVB[[3]]$probLeft, SVB[[4]]$probLeft,
                                           UVB[[1]]$probLeft, UVB[[2]]$probLeft, UVB[[3]]$probLeft, UVB[[4]]$probLeft,
                                           MCMC[[1]]$probLeft, MCMC[[2]]$probLeft),
                              probRight = c(SVB[[1]]$probRight, SVB[[2]]$probRight, SVB[[3]]$probRight, SVB[[4]]$probRight,
                                           UVB[[1]]$probRight, UVB[[2]]$probRight, UVB[[3]]$probRight, UVB[[4]]$probRight,
                                           MCMC[[1]]$probRight, MCMC[[2]]$probRight),
                              h = rep(1:H, 10),
                              startLane = data[sSeq[s], 6],
                              endLane = data[sSeq[s]+1:H, 6],
                              model = rep(c('IH', 'IH', 'CH', 'CH', 'IH', 'IH', 'CH', 'CH', 'IH', 'CH'), rep(H, 10)),
                              inference = rep(c(rep('SVB', 4), rep('UVB', 4), rep('MCMC', 2)), rep(H, 10)),
                              S = sSeq[s],
                              id = id)
                   )
  # Homogenous model results
  homogResults <-  MCMCDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                            N = 200,
                            H = H,
                            grid = grid,
                            MCMCdraws = homogDraws)
  
  results <- rbind(results,
                   data.frame(logscore = homogResults$logscore,
                              xCDF = homogResults$xCDF,
                              dist = homogResults$mapDist,
                              probLeft = homogResults$probLeft,
                              probRight = homogResults$probRight,
                              h = 1:H,
                              startLane = data[sSeq[s], 6],
                              endLane = data[sSeq[s]+1:H, 6],
                              model = 'Homogenous',
                              inference = 'MCMC',
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
  v0 <- rep(data[sSeq[s] + 1, 3], 9) 
  
  laneSeq <- seq(-5.5 * laneSize, 5.5 * laneSize, laneSize)
  currentLane <- max(which(laneSeq < data[sSeq[s], 4]))
  
  for(h in 1:H){
    v0 <- v0 + c(rep(aConst, 3), rep(aAvg, 3), rep(0, 3))
    x0 <- x0 + v0 * rep(c(cos(dConst + pi/2), cos(pi/2), cos(pi/2 + dAvg)), 3)
    y0 <- y0 + v0 * rep(c(sin(dConst + pi/2), sin(pi/2), sin(pi/2 + dAvg)), 3)
    dist <- sqrt((x0 - sum(data[sSeq[s] + 1:h, 4]))^2 + (y0 - sum(data[sSeq[s] + 1:h, 5]))^2)
    results <- rbind(results,
                     data.frame(logscore = NA,
                                xCDF = NA,
                                dist = dist,
                                probLeft = x0 < laneSeq[currentLane],
                                probRight = x0 > laneSeq[currnetLane + 1],
                                h = h,
                                startLane = data[sSeq[s], 6],
                                endLane = data[sSeq[s]+h, 6],
                                model = paste('Naive', 1:9),
                                inference = NA,
                                S = sSeq[s],
                                id = id))
  }
  
  print(paste(i, s))
}

write.csv(results, paste0('eval/car', id, '.csv'), row.names=FALSE)
