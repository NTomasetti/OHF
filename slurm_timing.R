rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(repenv)
set.seed(1000 + i)

library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
library(RcppEigen, lib.loc = 'packages')
library(rstan, lib.loc = 'packages')
source('slurm_RFuns.R')
sourceCpp('slurm_cppFuns.cpp')

ids <- readRDS('idTest.RDS')
datatest <- readRDS('dataTest.RDS')

data <- datatest[[i]]
id <- ids[i]

CHprior <- readRDS('CHFit.RDS')




minT <- 50
maxT <- 400

results <- data.frame()
lagsA <- 4
lagsD <- 4
K <- mix <- 6
dim <- 2 + lagsA + lagsD

# VB Prior Distributions: 1) IH, 2) CH
prior <- list()
prior[[1]] <- c(-5, -5, rep(0, 8), c(chol(diag(10, 10))))
prior[[2]] <- CHprior

starting <- list(matrix(c(prior[[2]]$mean[1:10], diag(0.5, 10)), ncol = 1),
                 matrix(c(prior[[2]]$mean, rep(log(0.5), 10*K), rep(1, K)), ncol = 1))


fitInitial <- fitVB(data[1:minT, 1:2], prior, starting, dim = 10, mix = K, time = TRUE)

S <- c(1, 5, 10, 20)

results <- data.frame(timing = fitInitial[[2]],
                      model = c('IH', 'IH', 'CH', 'CH'),
                      inference =  c('SVB-Single', 'SVB-Mixture', 'SVB-Single', 'SVB-Mixture', 'UVB-Single', 
                                     'UVB-Mixture', 'UVB-Single', 'UVB-Mixture'),
                      k = rep(S, rep(8, 4)),
                      S = minT,
                      id = id)

for(k in 1:4){
  
  sSeq <- seq(minT+S[k], maxT, S[k])
  
  # Incrementally add data to VB fits
  for(s in seq_along(sSeq)){
    if(sSeq[s] > nrow(data)){
      break
    }
    
    # Offline (Standard VB) model fits
    fitOffline <- fitVB(data[1:sSeq[s], 1:2], prior, starting, dim = 10, mix = K, time = TRUE)
    
    if(s == 1){
      fitOnline <- fitUVB(data[(minT-(max(lagsA, lagsD) - 1)):sSeq[s], 1:2], fitInitial[[1]], starting, dim = 10, mix = K, time = TRUE)
    } else {
      fitOnline <- fitUVB(data[(sSeq[s-1]-(max(lagsA, lagsD) - 1)):sSeq[s], 1:2], fitOnline[[1]], starting, dim = 10, mix = K, time = TRUE)
    }
   
    
    # Grab logscores etc. for heterogenous models
    results <- rbind(results,
                     data.frame(timing = c(fitOffline[[2]], fitOnline[[2]]),
                                model = c('IH', 'IH', 'CH', 'CH'),
                                inference = c('SVB-Single', 'SVB-Mixture', 'SVB-Single', 'SVB-Mixture', 'UVB-Single', 
                                                  'UVB-Mixture', 'UVB-Single', 'UVB-Mixture'),
                                k = S[k],
                                S = sSeq[s],
                                id = id)
    )
  }
}

write.csv(results, paste0('timing/car', id, '.csv'), row.names=FALSE)

