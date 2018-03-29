library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
sourceCpp('MCMCFuns.cpp')

# An MCMC algorithn for the clustered hierarchical model with a K component gaussian prior distribution
# Each theta_i is drawn with metropolis hastings, 
# Prior components (beta) are drawn with gibbs from the posterior conditionals
mixtureMCMC <- function(data, reps, draw, hyper, K = 2, thinB = 10, thinT = 100, lagsA = 1, lagsD = 1, includeRho = FALSE, stepsize = 0.01){
  N <- length(data)
  stepsize <- rep(stepsize, N)
  accept <- rep(0, N)
  
  dim <- 2 + lagsA + lagsD + includeRho
  #set up storage for saved draws
  nSaveB <- floor(reps / thinB)
  nSaveT <- floor(reps / thinT)
  saveDraws <- list(list())
  for(i in 1:K){
    saveDraws[[1]][[i]] <- list(mean = matrix(0, nSaveB, dim), varInv = array(0, dim = c(dim, dim, nSaveB)))
  }
  for(i in 1:N){
    saveDraws[[i+1]] <- list(theta = matrix(0, nSaveT, dim), k = rep(0, nSaveT), pi = matrix(0, nSaveT, K))
  }
  # changing MH acceptance rate
  alpha <- - qnorm(0.234/2)
  stepsizeCons <- (1 - 1/dim) * sqrt(2*pi) * exp(alpha^2/2) / (2 * alpha)  +  1 / (dim * 0.234 * (1 - 0.234))
  
  for(i in 2:reps){
    # timing
    if(i == 50){
      startTime <- Sys.time()
    } else if(i == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
    }
    # theta_i
    for(j in 1:N){
      candidate <- draw[[j+1]]$theta +  stepsize[j] * rnorm(dim)
      group <- draw[[j+1]]$k
      canDens <- nlogDensity(data[[j]], candidate, draw[[1]][[group]]$mean, draw[[1]][[group]]$varInv, lagsA, lagsD, includeRho)
      oldDens <- nlogDensity(data[[j]], draw[[j+1]]$theta, draw[[1]][[group]]$mean, draw[[1]][[group]]$varInv, lagsA, lagsD, includeRho)
      
      ratio <- exp(canDens - oldDens)
      
      c <- stepsize[j] * stepsizeCons
      if(runif(1) < ratio){
        accept[j] <- accept[j] + 1
        draw[[j+1]]$theta <- candidate
        stepsize[j] <- stepsize[j] + c * (1 - 0.234) / (28 + i)
      } else {
        stepsize[j] <- stepsize[j] - c * 0.234 / (28 + i)
      }
    }
    # k_i
    for(j in 1:N){
      p <- numeric(K)
      for(k in 1:K){
        p[k] <-  log(draw[[j+1]]$pi[k])  +  0.5 * log(det(draw[[1]][[k]]$varInv)) -
          0.5 * (draw[[j+1]]$theta - draw[[1]][[k]]$mean) %*% draw[[1]][[k]]$varInv %*%
          (draw[[j+1]]$theta - draw[[1]][[k]]$mean)
      }
      p <- p - max(p) 
      p <- exp(p) / sum(exp(p))
      draw[[j+1]]$k <- base::sample(1:K, 1, prob=p)
    }
    #pi_i
    for(j in 1:N){
      group <- rep(0, K)
      group[draw[[j+1]]$k] <- 1
      draw[[j+1]]$pi <- c(rdirichlet(1, hyper$alpha + group))
    }

    # thetaHat_k
    sumK <- rep(0, K)
    sumTheta <- matrix(0, dim, K)
    for(j in 1:N){
      sumK[draw[[j+1]]$k] <- sumK[draw[[j+1]]$k] + 1
      sumTheta[,draw[[j+1]]$k] <-  sumTheta[,draw[[j+1]]$k] + draw[[j+1]]$theta
    }
    for(k in 1:K){
      var <- hyper[[k]]$varInv + sumK[k] * draw[[1]][[k]]$varInv
      var <- solve(var)
      mean <- var %*% hyper[[k]]$varInv %*% hyper[[k]]$mean  +  var %*% draw[[1]][[k]]$varInv %*%  sumTheta[,k]
      
      draw[[1]][[k]]$mean <- c(rmvnorm(1, mean, var))
      vardf <- hyper[[k]]$v 
      scaleMat <- hyper[[k]]$scale
      for(j in 1:N){
        if(draw[[j+1]]$k == k){
          vardf <- vardf + 1
          scaleMat <- scaleMat + outer(draw[[j+1]]$theta - draw[[1]][[k]]$mean, draw[[j+1]]$theta - draw[[1]][[k]]$mean)
        }
      }
      draw[[1]][[k]]$varInv <- rWishart(1, vardf, solve(scaleMat))[,,1]
    }
    # save draws
    if(i %% thinB == 0){
      for(k in 1:K){
        saveDraws[[1]][[k]]$mean[i/thinB,] <- draw[[1]][[k]]$mean
        saveDraws[[1]][[k]]$varInv[,,i/thinB] <- draw[[1]][[k]]$varInv
      }
    }
    if(i %% thinT == 0){
      for(j in 1:N){
        saveDraws[[j+1]]$theta[i/thinT, ] <- draw[[j+1]]$theta
        saveDraws[[j+1]]$k[i/thinT] <- draw[[j+1]]$k
        saveDraws[[j+1]]$pi[i/thinT, ] <- draw[[j+1]]$pi
      }
    }
    # print progress
    if(i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else if(mins > 5) {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      } else if(mins > 1){
        secs <- floor(mins * 60) %% 60
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', floor(mins), ' minutes and ', secs, ' seconds.'))
      } else {
        secs <- floor(mins * 60)
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', secs, ' seconds.'))
      }
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}

# An MCMC algorithm for the homogenous model
# Variance components are drawn with metropolis hastings while ar components are drawn from gibbs posterior conditionals
homogMCMC <- function(data, reps, draw, hyper, thin = 1, lagsA = 1, lagsD = 1, includeRho = FALSE, stepsize = 0.01){
  N <- length(data)
  # set up likelihood function and theta dimension
  dim <- 2 + lagsA + lagsD + includeRho
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  # Hyper Variance
  hyperVar <- diag(solve(hyper$varInv))

  # Various sums of data for the conditional distribution cross terms
  xSums <- matrix(0, 2 + lagsA + lagsD, 2 + lagsA + lagsD)
  lag <- max(lagsA, lagsD) + 1
  for(i in 1:N){
    T <- nrow(data[[i]])
    # A_t^2
    xSums[1, 1] <- xSums[1, 1] + sum(data[[i]][lag:T, 1]^2)
    # A_t * A_t-j
    for(j in 1:lagsA){
      xSums[1, 1+j] <- xSums[1, 1+j] + sum(data[[i]][lag:T, 1] * data[[i]][lag:T - j, 1])
    }
    # A_t * D_t
    xSums[1, 2+lagsA] <- xSums[1, 2+lagsA] + sum(data[[i]][lag:T, 1] * data[[i]][lag:T, 2])
    # A_t * D_t-j
    for(j in 1:lagsD){
      xSums[1, 2+lagsA+j] <- xSums[1, 2+lagsA+j] + sum(data[[i]][lag:T, 1] * data[[i]][lag:T - j, 2])
    }
    for(j in 1:lagsA){
      # A_t-j * A_t-k
      for(k in 1:j){
        xSums[1 + k, 1 + j] <- xSums[1+k, 1 + j] + sum(data[[i]][lag:T - j, 1] * data[[i]][lag:T - k, 1])
      }
      # A_t-j * D_t
      xSums[1+j, 2+lagsA] <- xSums[1+j, 2+lagsA] + sum(data[[i]][lag:T -j, 1] * data[[i]][lag:T, 2])
      # A_t-j * D_t-k
      for(k in 1:lagsD){
        xSums[1 + j, 2 + lagsA + k] <- xSums[1+j, 2 + lagsA + k] + sum(data[[i]][lag:T - j, 1] * data[[i]][lag:T - k, 2])
      }
    }
    #D_t * D_t
    xSums[2 + lagsA, 2 + lagsA] <- xSums[2 + lagsA, 2 + lagsA] + sum(data[[i]][lag:T, 2]^2)
    #D_t * D_t-j
    for(j in 1:lagsD){
      xSums[2+lagsA, 2+lagsA+j] <- xSums[2+lagsA, 2+lagsA+j] + sum(data[[i]][lag:T, 2] * data[[i]][lag:T - j, 2])
    }
    #D_t-j * D_t-k
    for(j in 1:lagsD){
      # A_t-j * A_t-k
      for(k in 1:j){
        xSums[2 + lagsA + k, 2 + lagsA + j] <- xSums[2 + lagsA + k, 2 + lagsA + j] + sum(data[[i]][lag:T - j, 2] * data[[i]][lag:T - k, 2])
      }
    }
  }
  xSums[lower.tri(xSums)] <- t(xSums)[lower.tri(xSums)]
  
 
  # Metropolis Hastings variables
  MHvar <- ifelse(includeRho, 1:3, 1:2)
  
  for(i in 2:reps){
    # timing
    if(i == 50){
      startTime <- Sys.time()
    } else if(i == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    }
    # log_sigma_squared and rho
    candidate <- draw
    candidate[MHvar] <- candidate[MHvar] + stepsize * rnorm(length(MHvar))
    canDens <- 0
    oldDens <- 0
    for(j in 1:N){
      canDens <- canDens + nlogDensity(data[[j]], candidate, hyper$mean, hyper$varInv, lagsA, lagsD, includeRho)
      oldDens <- oldDens + nlogDensity(data[[j]], draw, hyper$mean, hyper$varInv, lagsA, lagsD, includeRho)
    }
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      draw <- candidate
      stepsize <- stepsize + c * (1 - 0.44) / (18 + i)
    } else {
      stepsize <- stepsize - c * 0.44 / (18 + i)
    }
    
    # Transform to variance params
    sigmaA <- sqrt(exp(draw[1]))
    sigmaD <- sqrt(exp(draw[2]))
    if(includeRho){
      rho <- 2 / (1 + exp(-draw[3])) - 1
    } else {
      rho <- 0
    }
    
    # Phi parameters
    for(k in 1:lagsA){
      # a_t * a_t terms
      meanNumer <- xSums[1, 1 + k] 
      for(l in 1:lagsA){
        if(l != k){
          meanNumer <- meanNumer - draw[2 + includeRho + l] * xSums[1 + k, 1 + l] 
        }
      }
      meanNumer <- meanNumer * sigmaD * hyperVar[2 + includeRho + k]
      # a_t * d_t terms
      meanNumer <- meanNumer - rho * sigmaA * hyperVar[2 + includeRho + k] * xSums[1 + k, 2 + lagsA]
      for(l in 1:lagsD){
        meanNumer <- meanNumer + rho * sigmaA * hyperVar[2 + includeRho + k] * draw[2 + includeRho + lagsA + l] * xSums[1 + k, 2 + lagsA + l]
      }
      meanDenom <- xSums[1+k, 1+k] * sigmaD * hyperVar[2 + includeRho + k] + sigmaA^2 * sigmaD * (1 - rho^2)
      varNumer <- (1 - rho^2) * sigmaA^2 * sigmaD * hyperVar[2 + includeRho + k]
      draw[2 + includeRho + k] <- rnorm(1, meanNumer / meanDenom, sqrt(varNumer / meanDenom))
    }
    
    # Gamma parameters
    for(k in 1:lagsD){
      # a_t * a_t terms
      meanNumer <- xSums[lagsA + 2, 2 + lagsA + k] 
      for(l in 1:lagsD){
        if(l != k){
          meanNumer <- meanNumer - draw[2 + includeRho + lagsA + l] * xSums[lagsA + 2 + k, lagsA + 2 + l] 
        }
      }
      meanNumer <- meanNumer * sigmaA * hyperVar[2 + includeRho + lagsA +  k]
      # a_t * d_t terms
      meanNumer <- meanNumer - rho * sigmaD * hyperVar[2 + includeRho + lagsA + k] * xSums[1, 2 + lagsA + k]
      for(l in 1:lagsA){
        meanNumer <- meanNumer + rho * sigmaD * hyperVar[2 + includeRho + lagsA + k] * draw[2 + includeRho + l] * xSums[2 + lagsA + k, 1 + l]
      }
      meanDenom <- xSums[2 + lagsA + k, 2 + lagsA + k] * sigmaA * hyperVar[2 + includeRho + lagsA + k] + sigmaD^2 * sigmaA * (1 - rho^2)
      varNumer <- (1 - rho^2) * sigmaD^2 * sigmaA * hyperVar[2 + includeRho + lagsA + k]
      draw[2 + includeRho + lagsA + k] <- rnorm(1, meanNumer / meanDenom, sqrt(varNumer / meanDenom))
    }
    
    # save draws
    if(i %% thin == 0){
      saveDraws[i/thin,] <- draw
    }
    # print progress
    if(i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else if(mins > 5) {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      } else if(mins > 1){
        secs <- floor(mins * 60) %% 60
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', floor(mins), ' minutes and ', secs, ' seconds.'))
      } else {
        secs <- floor(mins * 60)
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', secs, ' seconds.'))
      }
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}

# An MCMC algorithm for the independent model
# The entire theta is drawn at the same time with metropolis hastings
# The model can be a var structure or a standard ar(lags).
indepMCMC <- function(data, reps, draw, hyper, thin = 1, stepsize = 0.01, lagsA = 1, lagsD = 1, VAR = FALSE){
  # set up likelihood function and theta dimension
  if(VAR){
    dim <- 3 + 2 * lagsA + 2 * lagsD
    likelihood <- nlogDensityVAR
  } else {
    dim <- 3 + lagsA + lagsD
    likelihood <- nlogDensity
  }
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  oldDens <- likelihood(data, draw, hyper$mean, hyper$varInv, lagsA, lagsD)
  
  for(i in 2:reps){
    candidate <- draw
    candidate <- candidate + stepsize * rnorm(dim)
    canDens <- likelihood(data, candidate, hyper$mean, hyper$varInv, lagsA, lagsD)
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      draw <- candidate
      stepsize <- stepsize + c * (1 - 0.44) / (18 + i)
      oldDens <- canDens
    } else {
      stepsize <- stepsize - c * 0.44 / (18 + i)
    }
    
    # save draws
    if(i %% thin == 0){
      saveDraws[i/thin,] <- draw
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}
