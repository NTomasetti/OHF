
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


singleMCMCallMH <- function(data, reps, draw, hyper, lagsA = 4, lagsD = 4, thin = 1, stepsize = 0.01, mix = FALSE){
  # set up theta dimension
  dim <- 2 + lagsA + lagsD
  # count accepted draws
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  if(mix){
    oldDens <- nMixLogDens(data, draw, hyper$mean, hyper$varInv, hyper$weights, lagsA, lagsD)
  } else {
    oldDens <- nlogDensity(data, draw, hyper$mean, hyper$varInv, lagsA, lagsD)
  }
  for(i in 2:reps){
    
    
    stat <- FALSE
    while(!stat){
      candidate <- draw + stepsize * rnorm(dim)
      poly1 <- polyroot(c(1, -candidate[3:(2 + lagsA)]))
      poly2 <- polyroot(c(1, -candidate[(2 + lagsA + 1):(dim)]))
      if(all(Mod(poly1) > 1) & all(Mod(poly2) > 1)){
        stat <- TRUE
      }
    }
    
    
    if(mix){
      canDens <- nMixLogDens(data, candidate, hyper$mean, hyper$varInv, hyper$weights, lagsA, lagsD)
    } else {
      canDens <- nlogDensity(data, candidate, hyper$mean, hyper$varInv, lagsA, lagsD)
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

carsVBScore <- function(data, lambda, model, dimTheta, K, S = 50, maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01, ...){
  
  dimLambda <- nrow(lambda)
  sobol <- sobol_points(100+6*S, dimTheta)
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
    z <- lambda[2*dimTheta*K + 1:K]
    pi <- exp(z) / sum(exp(z))
    unif <- shuffle(sobol)
    s <- 0
    try <- 0
    epsilon <- qnorm(unif[101:(100+6*S), ])
    epsilon[epsilon < -3] = -3
    epsilon[epsilon > 3] = 3
    while(s < S){
      k <- sample(1:K, 1, prob=pi)
      Qmean <- lambda[(k-1)*dim + 1:dim]
      Qsd <- exp(lambda[dim*mix + (k-1)*dim + 1:dim])
      theta <- c(Qmean + Qsd * epsilon[try+1,])
      derivs <- model(data, lambda, theta, mix = K, ...)
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

fitVB <- function(data, prior, starting, dim, mix){
  # Fit Standard VB (ie. Offline)
  fit <- list()
  
  # IH / Single Component Approx
  mean <- prior[[1]][1:dim]
  u <- matrix(prior[[1]][dim + 1:dim^2], dim)
  linv <- solve(t(u))
  fit[[1]] <- carsVB(data, starting[[1]], model = singlePriorSingleApprox, S = 10, dimTheta = dim, mean = mean, Linv = linv,  lagsA = 4, lagsD =4)$lambda
  
  # IH / Mixture
  fit[[2]] <- carsVBScore(data, starting[[2]], model = singlePriorMixApprox, dimTheta = dim, K = mix, mean = mean, Linv = linv, lagsA = 4, lagsD =4)$lambda
  
  # CH / Single Component Approx
  fit[[3]] <- carsVB(data, starting[[1]], model = mixPriorSingleApprox, S = 10, dimTheta = dim, mean = prior[[2]]$mean, Linv = prior[[2]]$linv,
                     dets = prior[[2]]$dets, weights = prior[[2]]$weights, lagsA = 4, lagsD =4)$lambda
  
  # CH / Mixture
  fit[[4]] <- carsVBScore(data, starting[[2]], model = mixPriorMixApprox, dimTheta = dim, K = mix, mean = prior[[2]]$mean, SigInv = prior[[2]]$varInv,
                          dets = prior[[2]]$dets, weights = prior[[2]]$weights, lagsA = 4, lagsD = 4)$lambda
  
  return(fit)
}

fitUVB <- function(data, prior, starting, dim, mix){
  fit <- list()
  
  # IH / Single Component Approx
  mean <- prior[[1]][1:dim]
  u <- matrix(prior[[1]][dim + 1:dim^2], dim)
  linv <- solve(t(u))
  fit[[1]] <- carsVB(data, starting[[1]], model = singlePriorSingleApprox, dimTheta = dim, mean = mean, Linv = linv, lagsA = 4, lagsD = 4)$lambda
  
  # IH / Mixture
  mean <- matrix(prior[[2]][1:(dim*mix)], dim)
  siginv <- array(0, dim = c(dim, dim, mix))
  dets <- NULL
  for(k in 1:mix){
    sd <- exp(prior[[2]][dim*mix + (k-1)*dim + 1:dim])
    var <- diag(sd^2)
    siginv[,,k] <- solve(var)
    dets <- c(dets, 1 / prod(sd))
  }
  weights <- prior[[2]][2*dim*mix + 1:mix]
  weights <- exp(weights) / sum(exp(weights))
  fit[[2]] <- carsVBScore(data, starting[[2]], model = mixPriorMixApprox, dimTheta = dim, K = mix, mean = mean, SigInv = siginv, 
                          dets = dets, weights = weights, lagsA = 4, lagsD = 4)$lambda
  
  # CH / Single Component Approx
  mean <- prior[[3]][1:dim]
  u <- matrix(prior[[3]][dim + 1:dim^2], dim)
  linv <- solve(t(u))
  fit[[3]] <- carsVB(data, starting[[1]], model = singlePriorSingleApprox, dimTheta = dim, mean = mean, Linv = linv, lagsA = 4, lagsD = 4)$lambda
  
  # CH / Mixture
  mean <- matrix(prior[[4]][1:(dim*mix)], dim)
  siginv <- array(0, dim = c(dim, dim, mix))
  dets <- NULL
  for(k in 1:mix){
    sd <- exp(prior[[4]][dim*mix + (k-1)*dim + 1:dim])
    var <- diag(sd^2)
    siginv[,,k] <- solve(var)
    dets <- c(dets, 1 / prod(sd))
  }
  weights <- prior[[4]][2*dim*mix + 1:mix]
  weights <- exp(weights) / sum(exp(weights))
  fit[[4]] <- carsVBScore(data, starting[[2]], model = mixPriorMixApprox, dimTheta = dim, K = mix, mean = mean, SigInv = siginv, 
                          dets = dets, weights = weights, lagsA = 4, lagsD = 4)$lambda
  
  return(fit)
}
  
autocov <- function(phi, maxLag, sigma2){
  
  p <- length(phi)
  
  maxLagp1 <- maxLag + 1
  r <- p + 1
  b <- numeric(r)
  phi2 <- numeric(3 * r)
  phi2[r] <- -1
  phi2[r + 1:p] <- phi
  
  b[1] <- b[1] + 1
  
  a <- matrix(numeric(r^2), ncol = r)
  for (i in 1:r){
    for (j in 1:r) {
      if (j == 1) 
        a[i, j] <- phi2[r + i - 1]
      else 
        a[i, j] <- phi2[r + i - j] + phi2[r + i + j - 2]
    }
  } 
  g <- solve(a, -b)
  if (length(g) <= maxLag) {
    g <- c(g, numeric(maxLag - r))
    for (i in (r + 1):maxLagp1) {
      g[i] <- phi %*% g[i - 1:p]
    }
  }
  
  sigma2 * g[1:maxLagp1]
}

MCMCDens <- function(data, H, grid, MCMCdraws, lagsA = 1, lagsD = 1){
  
  samples <- nrow(MCMCdraws)
  # Calcualte prediction variance using ltsa
  predVar <- array(0, dim = c(samples, H, 2))
  for(i in 1:samples){
    covA <- autocov(phi = c(MCMCdraws[i, 3:6]), maxLag = H, sigma2 = exp(MCMCdraws[i, 1]))
    predVar[i, , 1] <- ltsa::PredictionVariance(covA, maxLead = H)
    covD <- autocov(phi = c(MCMCdraws[i, 7:10]), maxLag = H, sigma2 = exp(MCMCdraws[i, 2]))
    predVar[i, , 2] <- ltsa::PredictionVariance(covD, maxLead = H)
  }
  
  # Evaluate forecast density up to H steps ahead along the grid
  adDens <-  evalDensity(data, MCMCdraws, nrow(grid), H, grid, predVar, lagsA, lagsD)
  
  # Set up results
  results <- data.frame()
  adGrid <- expand.grid(a = grid[,1], d = grid[,2])
  starting <- max(lagsA, lagsD)
  mapX <- 0
  mapY <- 0
  for(h in 1:H){
    # Set up grid of all combinations of a and d
    fullGrid <- cbind(adGrid, expand.grid(aD = adDens[,1,h], dD = adDens[,2,h]))
    fullGrid$v <- fullGrid$a + data[starting + h, 3]
    # Transform a / d density to x / y density
    fullGrid$dens <- fullGrid$aD * fullGrid$dD / fullGrid$v
    fullGrid$x <- fullGrid$v * cos(fullGrid$d + 3.141598/2)
    fullGrid$y <- fullGrid$v * sin(fullGrid$d + 3.141598/2) 
    # Calculate distance to true value
    fullGrid$dist <- sqrt((fullGrid$x - data[starting + h, 4] + data[starting + h - 1, 4])^2 + (fullGrid$y - data[starting + h, 5] + data[starting + h - 1, 5])^2)
    
    # Logscore = Log-density at location that minimises distance
    logscore <- log(fullGrid$dens[which.min(fullGrid$dist)])
    # Map = Location that maximises density
    mapX <- mapX + fullGrid$x[which.max(fullGrid$dens)]
    mapY <- mapY + fullGrid$y[which.max(fullGrid$dens)]
    # Distance from (cumulative) map to actual location
    mapDist <- sqrt((mapX - data[starting+h, 4] + data[starting, 4])^2 + (mapY - data[starting + h, 5] + data[starting, 5])^2)
    
    results <- rbind(results, data.frame(h = h, logscore = logscore, mapDist = mapDist))
    
  }
  results
}

rtruncnorm <- function(n, mu, sd, min, max){
  if(min >= max){
    stop('min must be smaller than max')
  }
  minP <- pnorm(min, mu, sd)
  maxP <- pnorm(max, mu, sd)
  if(round(minP, 6) == 1 | round(maxP, 6) == 0){
    runif(n, min, max)
  } else {
    u <- runif(n, minP, maxP)
    qnorm(u, mu, sd)
  }

}

VBDens <- function(data, fit, grid, H, dim, K, isMix, lagsA = 1, lagsD = 1, samples = 500){
  # This works the same as the MCMC density evaluation
  # But first we need to sample from the VB density to feed this into evalDensity()
  N <- nrow(grid)
  draws <- matrix(0, samples, 10)
  if(!isMix){
    mean <- fit[1:dim]
    L <- t(matrix(fit[dim + 1:dim^2], dim))
    weights <- 1
    
    for(i in 1:samples){
      
      stat <- FALSE
      try <- 0
      while(!stat){
        try <- try + 1
        candidate <- mean + L %*% rnorm(10)
        poly1 <- polyroot(c(1, -candidate[3:(2 + lagsA)]))
        poly2 <- polyroot(c(1, -candidate[(2 + lagsA + 1):(dim)]))
        if(all(Mod(poly1) > 1) & all(Mod(poly2) > 1)){
          stat <- TRUE
        }
        if(!stat & try > 5){
          Sigma <- L %*% t(L)
          if(!all(Mod(poly1) > 1)){
            if(lagsA == 1){
              candidate[3] <- rtruncnorm(1, mean[3], sqrt(Sigma[3, 3]), -1, 1)
            } else {
              candidate[3] <- rtruncnorm(1, mean[3], sqrt(Sigma[3, 3]), -2, 2)
              candidate[4] <- rtruncnorm(1, mean[4], sqrt(Sigma[4, 4]), -1, 1 - abs(candidate[3]))
              if(lagsA > 2){
                candidate[5:(2 + lagsA)] <- 0
              }
            }
          }
          if(!all(Mod(poly2) > 1)){
            if(lagsD == 1){
              candidate[3 + lagsA] <- rtruncnorm(1, mean[3 + lagsA], sqrt(Sigma[3 + lagsA, 3 + lagsA]), -1, 1)
            } else {
              candidate[3 + lagsA] <- rtruncnorm(1, mean[3 + lagsA], sqrt(Sigma[3 + lagsA, 3 + lagsA]), -2, 2)
              candidate[4 + lagsA] <- rtruncnorm(1, mean[4 + lagsA], sqrt(Sigma[4 + lagsA, 4 + lagsA]), -1, 1 - abs(candidate[3 + lagsA]))
              if(lagsD > 2){
                candidate[(5 + lagsA):dim] <- 0
              }
            }
          }
          stat <- TRUE
        }
      }
      draws[i, ] <- candidate
    }
    
  } else {
    z <- fit[2*dim*K + 1:K]
    pi <- exp(z) / sum(exp(z))
    weights <- cumsum(pi)
    mean <- matrix(fit[1:(dim*K)], ncol = K)
    L <- array(0, dim = c(dim, dim, K))
    for(k in 1:K){
      sd <- exp(fit[dim*K + (k-1)*dim + 1:dim])
      L[,,k] <- diag(sd)
    }
    for(i in 1:samples){
      group <- sample(1:K, 1, prob = pi)
      stat <- FALSE
      try <- 0
      while(!stat){
        try <- try + 1
        candidate <- mean[,group] + L[,,group] %*% rnorm(10)
        poly1 <- polyroot(c(1, -candidate[3:(2 + lagsA)]))
        poly2 <- polyroot(c(1, -candidate[(2 + lagsA + 1):(dim)]))
        if(all(Mod(poly1) > 1) & all(Mod(poly2) > 1)){
          stat <- TRUE
        }
        if(!stat & try > 5){
          Sigma <- L[,,group] %*% t(L[,,group])
          if(!all(Mod(poly1) > 1)){
            if(lagsA == 1){
              candidate[3] <- rtruncnorm(1, mean[3, group], sqrt(Sigma[3, 3]), -1, 1)
            } else {
              candidate[3] <- rtruncnorm(1, mean[3, group], sqrt(Sigma[3, 3]), -2, 2)
              candidate[4] <- rtruncnorm(1, mean[4, group], sqrt(Sigma[4, 4]), -1, 1 - abs(candidate[3]))
              if(lagsA > 2){
                candidate[5:(2 + lagsA)] <- 0
              }
            }
          }
          if(!all(Mod(poly2) > 1)){
            if(lagsD == 1){
              candidate[3 + lagsA] <- rtruncnorm(1, mean[3 + lagsA, group], sqrt(Sigma[3 + lagsA, 3 + lagsA]), -1, 1)
            } else {
              candidate[3 + lagsA] <- rtruncnorm(1, mean[3 + lagsA, group], sqrt(Sigma[3 + lagsA, 3 + lagsA]), -2, 2)
              candidate[4 + lagsA] <- rtruncnorm(1, mean[4 + lagsA, group], sqrt(Sigma[4 + lagsA, 4 + lagsA]), -1, 1 - abs(candidate[3 + lagsA]))
              if(lagsD > 2){
                candidate[(5 + lagsA):dim] <- 0
              }
            }
          }
          stat <- TRUE
        }
      }
      draws[i, ] <- candidate
    }
  }
  
  # From here this is the same as the MCMC version
  starting <- max(lagsA, lagsD)

  predVar <- array(0, dim = c(samples, H, 2))
  for(i in 1:samples){
    covA <- autocov(phi = c(draws[i, 3:6]), maxLag = H, sigma2 = exp(draws[i, 1]))
    predVar[i, , 1] <- ltsa::PredictionVariance(covA, maxLead = H)
    covD <- autocov(phi = c(draws[i, 7:10]), maxLag = H, sigma2 = exp(draws[i, 2]))
    predVar[i, , 2] <- ltsa::PredictionVariance(covD, maxLead = H)
  }
  
  adDens <- evalDensity(data, draws, nrow(grid), H, grid, predVar, lagsA, lagsD)
  results <- data.frame()
  adGrid <- expand.grid(a = grid[,1], d = grid[,2])
  mapX <- 0
  mapY <- 0
  for(h in 1:H){
    fullGrid <- cbind(adGrid, expand.grid(aD = adDens[,1,h], dD = adDens[,2,h]))
    fullGrid$v <- fullGrid$a + data[starting + h, 3]
    fullGrid$dens <- fullGrid$aD * fullGrid$dD / fullGrid$v
    fullGrid$x <- fullGrid$v * cos(fullGrid$d + 3.141598/2)
    fullGrid$y <- fullGrid$v * sin(fullGrid$d + 3.141598/2) 
    fullGrid$dist <- sqrt((fullGrid$x - data[starting + h, 4] + data[starting + h - 1, 4])^2 + (fullGrid$y - data[starting + h, 5] + data[starting + h - 1, 5])^2)
    
    logscore <- log(fullGrid$dens[which.min(fullGrid$dist)])
    mapX <- mapX + fullGrid$x[which.max(fullGrid$dens)]
    mapY <- mapY + fullGrid$y[which.max(fullGrid$dens)]
    mapDist <- sqrt((mapX - data[starting+h, 4] + data[starting, 4])^2 + (mapY - data[starting + h, 5] + data[starting, 5])^2)
	
    results <- rbind(results, data.frame(h = h, logscore = logscore, mapDist = mapDist))
    
  }
  results
}




