library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('MCMCFuns.cpp')

mixtureMCMC <- function(data, reps, draw, hyper, thin = 1, K = 2, error = 'gaussian', stepsize = 0.01){
  N <- length(data)
  stepsize <- rep(stepsize, N)
  accept <- rep(0, N)
  
  
  #set up likelihood function and theta dimension
  if(error == 'gaussian'){
    likelihood <- nlogDensity
    dim <- 6
  } else if(error == 't') {
    likelihood <- tlogDensity
    dim <- 8
  } else {
    stop('error must be gaussian or t')
  }
  #set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- list(list())
  for(i in 1:K){
    saveDraws[[1]][[i]] <- list(mean = matrix(0, nSave, dim), varInv = array(0, dim = c(dim, dim, nSave)))
  }
  for(i in 1:N){
    saveDraws[[i+1]] <- list(theta = matrix(0, nSave, dim), k = rep(0, nSave), pi = matrix(0, nSave, K))
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
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    }
    # theta_i
    for(j in 1:N){
      candidate <- draw[[j+1]]$theta +  stepsize[j] * rnorm(dim)
      group <- draw[[j+1]]$k
      canDens <- likelihood(data[[j]], candidate, draw[[1]][[group]]$mean, draw[[1]][[group]]$varInv)
      oldDens <- likelihood(data[[j]], draw[[j+1]]$theta, draw[[1]][[group]]$mean, draw[[1]][[group]]$varInv)
      
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
      draw[[j+1]]$pi <- c(MCMCpack::rdirichlet(1, hyper$alpha + group))
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
    if(i %% thin == 0){
      for(k in 1:K){
        saveDraws[[1]][[k]]$mean[i/thin,] <- draw[[1]][[k]]$mean
        saveDraws[[1]][[k]]$varInv[,,i/thin] <- draw[[1]][[k]]$varInv
      }
      for(j in 1:N){
        saveDraws[[j+1]]$theta[i/thin, ] <- draw[[j+1]]$theta
        saveDraws[[j+1]]$k[i/thin] <- draw[[j+1]]$k
        saveDraws[[j+1]]$pi[i/thin, ] <- draw[[j+1]]$pi
      }
    }
    # print progress
    if(i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      }
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}

homogMCMC <- function(data, reps, draw, hyper, thin = 1, error = 'gaussian', stepsize = 0.01){
  N <- length(data)
  # set up likelihood function and theta dimension
  if(error == 'gaussian'){
    likelihood <- nlogDensity
    dim <- 6
  } else if(error == 't') {
    likelihood <- tlogDensity
    dim <- 8
  } else {
    stop('error must be gaussian or t')
  }
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  
  # sums of data
  a11 <- a22 <- a12 <- a01 <- a02 <- 0
  d11 <- d22 <- d12 <- d01 <- d02 <- 0
  for(i in 1:N){
    T <- nrow(data[[i]])
    for(t in 3:T){
      a11 <- a11 + data[[i]][t-1, 1]^2
      a22 <- a22 + data[[i]][t-2, 1]^2
      a12 <- a12 + data[[i]][t-1, 1] * data[[i]][t-2, 1]
      a01 <- a01 + data[[i]][t-1, 1] * data[[i]][t, 1]
      a02 <- a02 + data[[i]][t-2, 1] * data[[i]][t, 1]
      d11 <- d11 + data[[i]][t-1, 2]^2
      d22 <- d22 + data[[i]][t-2, 2]^2
      d12 <- d12 + data[[i]][t-1, 2] * data[[i]][t-2, 2]
      d01 <- d01 + data[[i]][t-1, 2] * data[[i]][t, 2]
      d02 <- d02 + data[[i]][t-2, 2] * data[[i]][t, 2]
    }
  }
  # Metropolis Hastings variables
  MHvar <- 1:2
  if(error == 't'){
    MHvar <- c(MHvar, 7:8)
  }
  
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
    # log_sigma_squared and nu (if error ~ t)
    candidate <- draw
    candidate[MHvar] <- candidate[MHvar] + stepsize * rnorm(dim-4)
    canDens <- 0
    oldDens <- 0
    for(j in 1:N){
      canDens <- canDens + likelihood(data[[j]], candidate, hyper$mean, hyper$varInv)
      oldDens <- oldDens + likelihood(data[[j]], draw, hyper$mean, hyper$varInv)
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
    # phi 1
    meanNumer <- hyper$var[3] * (a01 - draw[4] * a12)  +  exp(draw[1]) * hyper$mean[3]
    meanDenom <- hyper$var[3] * a11  +  exp(draw[1])
    var <- exp(draw[1]) * hyper$var[3] / meanDenom
    draw[3] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    
    # phi 2
    meanNumer <- hyper$var[4] * (a02 - draw[3] * a12)  +  exp(draw[1]) * hyper$mean[4]
    meanDenom <- hyper$var[4] * a22  +  exp(draw[1])
    var <- exp(draw[1]) * hyper$var[4] / meanDenom
    draw[4] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    
    # gamma 1
    meanNumer <- hyper$var[5] * (d01 - draw[6] * d12)  +  exp(draw[2]) * hyper$mean[5]
    meanDenom <- hyper$var[5] * d11  +  exp(draw[2])
    var <- exp(draw[2]) * hyper$var[5] / meanDenom
    draw[5] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    
    # gamma 2
    meanNumer <- hyper$var[6] * (d02 - draw[5] * d12)  +  exp(draw[2]) * hyper$mean[6]
    meanDenom <- hyper$var[6] * d22  +  exp(draw[2])
    var <- exp(draw[2]) * hyper$var[6] / meanDenom
    draw[6] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    
    # save draws
    if(i %% thin == 0){
      saveDraws[i/thin,] <- draw
    }
    # print progress
    if(i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      }
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}

indepMCMC <- function(data, reps, draw, hyper, thin = 1, stepsize = 0.01, lags = 1){
  # set up likelihood function and theta dimension
  dim <- 3 + 2 * lags
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  
  for(i in 2:reps){
    candidate <- draw
    candidate <- candidate + stepsize * rnorm(dim)
    canDens <- nlogDensity(data, candidate, hyper$mean, hyper$varInv, lags)
    oldDens <- nlogDensity(data, draw, hyper$mean, hyper$varInv, lags)
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      draw <- candidate
      stepsize <- stepsize + c * (1 - 0.44) / (18 + i)
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

indepMCMCVAR <- function(data, reps, draw, hyper, thin = 1, stepsize = 0.01, lags = 1){
  # set up likelihood function and theta dimension
  dim <- 3 + 4 * lags
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  
  for(i in 2:reps){
    candidate <- draw
    candidate <- candidate + stepsize * rnorm(dim)
    canDens <- nlogDensityVAR(data, candidate, hyper$mean, hyper$varInv, lags)
    oldDens <- nlogDensityVAR(data, draw, hyper$mean, hyper$varInv, lags)
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      draw <- candidate
      stepsize <- stepsize + c * (1 - 0.44) / (18 + i)
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

indepMCMCexactVAR <- function(data, reps, draw, hyper, thin = 1, stepsize = 0.01, lags = 2){
  # set up likelihood function and theta dimension
  dim <- 3 + 4 * lags
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  
  T <- nrow(data)
  
  X <- matrix(0, T-lags, 2 + 4 * lags)
  X[,1] <- data[(lags+1):T, 1]
  X[,2*lags + 2] <- data[(lags+1):T, 2]
  
  for(i in 1:lags){
    X[,1 + i] <- data[(lags + 1 - i):(T - i), 1]
    X[,lags + 1 + i] <- sign(data[(lags + 1 - i):(T - i), 1]) * abs(data[(lags + 1 - i):(T-i), 2])
    X[,2*lags + 2 + i] <-  data[(lags + 1 - i):(T - i), 2]
    X[,3*lags + 2 + i] <- sign(data[(lags + 1 - i):(T - i), 2]) * abs(data[(lags + 1 - i):(T-i), 1])
  }
  
  Xsum <- matrix(0, 4*lags+2, 4*lags+2)
  for(i in 1:(4*lags + 2)){
    for(j in 1:(4*lags + 2)){
      Xsum[i, j] <- sum(X[,i] * X[,j])
    }
  }
  hyperVar <- diag(solve(hyper$varInv))
  
  for(iter in 2:reps){
    candidate <- draw
    candidate[1:3] <- candidate[1:3] + stepsize * rnorm(3)
    canDens <- nlogDensityVAR(data, candidate, hyper$mean, hyper$varInv, lags)
    oldDens <- nlogDensityVAR(data, draw, hyper$mean, hyper$varInv, lags)
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      draw <- candidate
      stepsize <- stepsize + c * (1 - 0.44) / (18 + i)
    } else {
      stepsize <- stepsize - c * 0.44 / (18 + i)
    }
    
    for(i in 1:(2 * lags)){
      # Draw Phis
      theta_num <- 3 + i
      meanNumer <- -Xsum[1, 1 + i]
      
      # mean at terms
      for(j in 1:(2 * lags)){
        if(j == i){
          next
        }
        meanNumer = meanNumer + draw[3 + j] * Xsum[1 + i, 1 + j]
      }
      meanNumer = meanNumer * hyperVar[theta_num] / exp(draw[1])
      # at dt cross product terms
      crossConst = 0.5 * hyperVar[theta_num] * (2 / (1 + exp(-draw[3])) - 1) / (sqrt(exp(draw[1])) * sqrt(exp(draw[2])))
      meanNumer = meanNumer + crossConst * Xsum[1 + i, 2 + 2 * lags]
      for(j in 1:(2 * lags)){
        meanNumer = meanNumer + crossConst * draw[3 + 2 * lags + j] * Xsum[1 + i, 2 + 2 * lags + j]
      }
      
      meanDenom <- hyperVar[theta_num] / exp(draw[1]) * Xsum[1 + i, 1 + i] + 1 + (1 - (2 / (1 + exp(-draw[3])) - 1)^2)
      var <- hyperVar[theta_num] * (1 + (1 - (2 / (1 + exp(-draw[3])) - 1)^2)) / meanDenom 
      
      draw[theta_num] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    }
    
    
    for(i in 1:(2 * lags)){
      # Draw Gammas
      theta_num <- 3 + 2 * lags + i
      meanNumer <- - Xsum[2 + 2 * lags, 2 + 2 * lags + i]
      
      # mean at terms
      for(j in 1:(2 * lags)){
        if(j == i){
          next
        }
        meanNumer = meanNumer + draw[3 + 2 * lags + j] * Xsum[2 + 2 * lags + i, 2 + 2 * lags + j]
      }
      meanNumer = meanNumer * hyperVar[theta_num] / exp(draw[2])
      # at dt cross product terms
      crossConst = 0.5 * hyperVar[theta_num] * (2 / (1 + exp(-draw[3])) - 1) / (sqrt(exp(draw[1])) * sqrt(exp(draw[2])))
      meanNumer = meanNumer + crossConst * Xsum[2 + 2 * lags + i, 1]
      for(j in 1:(2 * lags)){
        meanNumer = meanNumer + crossConst * draw[3 + j] * Xsum[2 + 2 * lags + i, 3 + j]
      }
      
      meanDenom <- hyperVar[theta_num] / exp(draw[2]) * Xsum[2 + 2 * lags + i, 2 + 2 * lags + i] + 1 + (1 - (2 / (1 + exp(-draw[3])) - 1)^2)
      var <- hyperVar[theta_num] * (1 + (1 - (2 / (1 + exp(-draw[3])) - 1)^2)) / meanDenom 
      
      draw[theta_num] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    }
    
    
    
    
    # save draws
    if(iter %% thin == 0){
      saveDraws[iter/thin,] <- draw
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}



