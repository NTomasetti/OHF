carsVB <- function(data, lambda, model, S = 25, maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01, 
                   dimTheta = 6,  ...){
  dimLambda <- length(lambda)
  sobol <- sobol_points(100+S, dimTheta)
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- numeric(dimLambda)
  V <- numeric(dimLambda)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    #if(any(is.na(lambda))){
    #  break
    #}
    grad <- matrix(0, dimLambda, S)
    eval <- numeric(S)
    q <- numeric(S)
    unif <- shuffle(sobol)
    epsilon <- qnorm(unif[101:(100+S), ])
    epsilon[epsilon < -3] = -3
    epsilon[epsilon > 3] = 3
    if(S == 1){
      logpj <- model(data, lambda, epsilon, ...)
      eval <- logpj$val
      grad <- logpj$grad
      q <- sum(dnorm(epsilon, log=TRUE))
      gradient <- grad
      gradientSq <- grad^2
      LB[iter] <- eval - q
    } else {
      for(s in 1:S){
        logpj <- model(data, lambda, epsilon[s,], ...)    
        eval[s] <- logpj$val
        grad[,s] <- logpj$grad
        q[s] <- sum(dnorm(epsilon[s,], log=TRUE))
      }
      eval[eval == -Inf] = NA
      gradient <- rowMeans(grad, na.rm = TRUE)
      gradientSq <- rowMeans(grad^2, na.rm=TRUE)
      LB[iter] <- mean(eval - q, na.rm=TRUE) 
    }
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    if(any(is.na(alpha * Mst / sqrt(Vst + e)))){
      print('Break')
      break
    }
    lambda <- lambda + alpha * Mst / sqrt(Vst + e)
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    iter <- iter + 1
  }
  return(list(lambda=lambda, LB = LB[1:min(iter-1, maxIter)], iter = min(maxIter, iter-1)))
}

singleMCMCallMH <- function(data, reps, draw, hyper, thin = 1, error = 'gaussian', stepsize = 0.01, mix = FALSE){
  # set up likelihood function and theta dimension
  dim <- 7
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  if(mix){
    oldDens <- nMixLogDens(data, draw, hyper$mean, hyper$varInv, hyper$weights)
  } else {
    oldDens <- nlogDensity(data, draw, hyper$mean, hyper$varInv)
  }
  for(i in 2:reps){
    candidate <- draw
    candidate <- candidate + stepsize * rnorm(dim)
    if(mix){
      canDens <- nMixLogDens(data, candidate, hyper$mean, hyper$varInv, hyper$weights)
    } else {
      canDens <- nlogDensity(data, candidate, hyper$mean, hyper$varInv)
    }
    
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      draw <- candidate
      oldDens <- canDens
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

carsVBScore <- function(data, lambda, model, S = 50, maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01, ...){
  
  dimLambda <- nrow(lambda)
  sobol <- sobol_points(100+6*S, 6)
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- rep(0, dimLambda)
  V <- rep(0, dimLambda)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold | iter < 100){
    if(iter > maxIter){
      break
    }
    grad <- matrix(0, dimLambda, S)
    eval <- numeric(S)
    z <- lambda[73:78]
    pi <- exp(z) / sum(exp(z))
    unif <- shuffle(sobol)
    s <- 0
    try <- 0
    epsilon <- qnorm(unif[101:(100+6*S), ])
    epsilon[epsilon < -3] = -3
    epsilon[epsilon > 3] = 3
    while(s < S){
      k <- sample(1:6, 1, prob=pi)
      Qmean <- lambda[(k-1)*6 + 1:6]
      Qsd <- exp(lambda[6*6 + (k-1)*6 + 1:6])
      theta <- c(Qmean + Qsd * epsilon[try+1,])
      derivs <- model(data, lambda, theta, ...)
      if(all(is.finite(derivs$grad)) & all(!is.na(derivs$grad)) & is.finite(derivs$val) & !is.na(derivs$val)){
        s <- s + 1
        eval[s] <- derivs$val
        grad[,s] <- derivs$grad
        if(s == S){
          gradient <- rowMeans(grad, na.rm = TRUE)
          gradientSq <- rowMeans(grad^2, na.rm = TRUE)
          LB[iter] <- mean(eval, na.rm = TRUE)
          break
        }
      }
      try <- try + 1
      if(try > 5*S){
        if(s > 1){
          gradient <- rowMeans(grad[,1:s], na.rm = TRUE)
          gradientSq <- rowMeans(grad[,1:s]^2, na.rm = TRUE)
          LB[iter] <- mean(eval[1:s], na.rm = TRUE)
        } else {
          LB[iter] <- LB[iter-1] - 1
        }
        break
      }
    }
    
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    if(any(is.na(alpha * Mst / sqrt(Vst + e)))){
      print('Break')
      break
    }
    lambda <- lambda + alpha * Mst / sqrt(Vst + e)
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    iter <- iter + 1
  }
  print(paste0('iter: ', min(iter-1, maxIter), ' ELBO: ', LB[min(iter-1, maxIter)]))
  return(list(lambda=lambda, LB = LB[1:min(iter-1, maxIter)], iter = min(maxIter, iter-1)))
}

fitCarMods <- function(data, prior, starting){
  # Fit Standard VB (ie. Offline)
  fit <- list()
  
  # IH / Single Component Approx
  mean <- prior[[1]][1:6]
  u <- matrix(prior[[1]][7:42], 6)
  linv <- solve(t(u))
  fit[[1]] <- carsVB(data, starting[[1]], lags = 2, model = singlePriorSingleApprox, mean = mean, Linv = linv, dimTheta = 6)$lambda
  
  #IH / Mixture
  fit[[2]] <- carsVBScore(data, starting[[2]], lags = 2, model = singlePriorMixApprox, mean = mean, Linv = linv)$lambda
  
  # CH / Single Component Approx
  mean <- matrix(prior[[2]][1:(6*6)], 6)
  siginv <- array(0, dim = c(6, 6, 6))
  dets <- NULL
  for(k in 1:6){
    sd <- exp(prior[[2]][6*6 + (k-1)*6 + 1:6])
    var <- diag(sd^2)
    siginv[,,k] <- solve(var)
    dets <- c(dets, 1 / prod(sd))
  }
  weights <- prior[[2]][73:78]
  weights <- exp(weights) / sum(exp(weights))
  start <- Sys.time()
  fit[[3]] <- carsVB(data, starting[[1]], lags = 2, model = mixPriorSingleApprox, mean = mean, sigInv = siginv, dets = dets, weights = weights, dimTheta = 6)$lambda
  
  # CH / Mixture
  fit[[4]] <- carsVBScore(data, starting[[2]], lags = 2, model = mixPriorMixApprox, mean = mean, Siginv = siginv, dets = dets, weights = weights)$lambda
  
  return(fit)
}
  
MCMCDens <- function(data, N, H, grid, MCMCdraws, laneSize){
  adDens <-  evalMCMCDensIndep(data, N, H, grid, MCMCdraws)

  results <- data.frame()
  adGrid <- expand.grid(a = grid[,1], d = grid[,2])

  mapX <- 0
  mapY <- 0
  laneSeq <- seq(-5.5 * laneSize, 5.5 * laneSize, laneSize)
  currentLane <- max(which(laneSeq < data[1, 4]))
  
  for(h in 1:H){
    fullGrid <- cbind(adGrid, expand.grid(aD = adDens[,1,h], dD = adDens[,2,h]))    
    fullGrid$v <- fullGrid$a + data[2 +	h, 3]
    fullGrid$dens <- fullGrid$aD * fullGrid$dD / fullGrid$v
    fullGrid$x <- fullGrid$v * cos(fullGrid$d + 3.141598/2)
    fullGrid$y <- fullGrid$v * sin(fullGrid$d + 3.141598/2)
    fullGrid$dist <- sqrt((fullGrid$x - data[2 + h, 4])^2 + (fullGrid$y - data[2 + h, 5])^2)


    cons <- sum(fullGrid$dens)
    xCDF <- sum(fullGrid$dens[fullGrid$x < data[2 + h, 4]]) / cons
    probLeft <- sum(fullGrid$dens[fullGrid$x < laneSeq[currentLane]]) / cons
    probRight <- sum(fullGrid$dens[fullGrid$x > laneSeq[currentLane + 1]]) / cons
    
    logscore <- log(fullGrid$dens[which.min(fullGrid$dist)])
    mapX <- mapX + fullGrid$x[which.max(fullGrid$dens)]
    mapY <- mapY + fullGrid$y[which.max(fullGrid$dens)]
    
    mapDist <- sqrt((mapX - sum(data[2+1:h, 4]))^2 + (mapY - sum(data[2 + 1:h, 5]))^2)
    
    
    results <- rbind(results, data.frame(h = h, xCDF = xCDF, logscore = logscore, mapDist = mapDist, probLeft = probLeft, probRight = probRight))

  }
  results
}

VBDens <- function(data, fit, grid, H, mix, laneSize){
  N <- nrow(grid)
  if(!mix){
    mean <- fit[1:6]
    L <- t(matrix(fit[7:42], 6))
    weights <- rep(1, 6)
  } else {
    z <- fit[73:78]
    pi <- exp(z) / sum(exp(z))
    weights <- cumsum(pi)
    mean <- fit[1:36]
    L <- matrix(0, 36, 6)
    for(k in 1:6){
      sd <- exp(fit[36 + (k-1)*6 + 1:6])
      L[1:6 + (k-1)*6, ] <- diag(sd)
    }
  }
  
  adDens <- evalVBDensIndep(data, mean, L, weights, N, H, grid, mix)
  results <- data.frame()
  adGrid <- expand.grid(a = grid[,1], d = grid[,2])
  
  mapX <- 0
  mapY <- 0
  currentLane <- max(which(laneSeq < data[1, 4]))
  for(h in 1:H){
    fullGrid <- cbind(adGrid, expand.grid(aD = adDens[,1,h], dD = adDens[,2,h]))
    fullGrid$v <- fullGrid$a + data[2 + h, 3]
    fullGrid$dens <- fullGrid$aD * fullGrid$dD / fullGrid$v
    fullGrid$x <- fullGrid$v * cos(fullGrid$d + 3.141598/2)
    fullGrid$y <- fullGrid$v * sin(fullGrid$d + 3.141598/2) 
    fullGrid$dist <- sqrt((fullGrid$x - data[2 + h, 4])^2 + (fullGrid$y - data[2 + h, 5])^2)    


    cons <- sum(fullGrid$dens)
    xCDF <- sum(fullGrid$dens[fullGrid$x < data[2 + h, 4]]) / cons
    probLeft <- sum(fullGrid$dens[fullGrid$x < laneSeq[currentLane]]) / cons
    probRight <- sum(fullGrid$dens[fullGrid$x > laneSeq[currentLane + 1]]) / cons
    
    logscore <- log(fullGrid$dens[which.min(fullGrid$dist)])
    
    mapX <- mapX + fullGrid$x[which.max(fullGrid$dens)]
    mapY <- mapY + fullGrid$y[which.max(fullGrid$dens)]
    
    mapDist <- sqrt((mapX - sum(data[2+1:h, 4]))^2 + (mapY - sum(data[2 + 1:h, 5]))^2)
	
    results <- rbind(results, data.frame(h = h, xCDF = xCDF, logscore = logscore, mapDist = mapDist, probLeft = probLeft, probRight = probRight))
    
  }
  results
}




