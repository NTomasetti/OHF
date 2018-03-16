library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
source('MCMCFuns.R')

read.csv('carsAug.csv') %>%
  select(ID, relX, relY, relA, relV, relD) %>%
  ungroup() -> carsAug
colnames(carsAug) <- c('ID', 'x', 'a', 'y', 'v', 'delta')

set.seed(1)
N <- 100
idSubset <- sample(unique(carsAug$ID), N)
data <- list()
for(i in 1:N){
  carsAug %>%
    ungroup() %>%
    filter(ID == idSubset[i]) %>%
    mutate(d = relD - pi/2) %>% 
    select(relA , d) %>%
    as.matrix() -> data[[i]]
}
saveRDS(data, 'MCMCData.RDS')
saveRDS(idSubset, 'trainingID.RDS')



HomogenoousModel {
  reps <- 50000
  
  #TODO: Clean up draws, hyper and posterior density plots for new model
  draws <- c(-5, -5, 0, 0, 0, 0)
  hyper <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(10, 6)))
  hyper$var <- diag(solve(hyper$varInv))
  homogDraws <- homogMCMC(data, reps, draws, hyper, thin = 10, stepsize = 0.01)
  saveRDS(homogDraws, 'homogMCMCN3000.RDS')
  
  homogDraws$draws %>%
    as.data.frame() %>%
    cbind(iter = 1:(reps/thin)) %>%
    mutate(V1 = exp(V1), V2 = exp(V2)) %>%
    rename(`sigma[epsilon]^{2}` = V1, `sigma[eta]^{2}` = V2, `phi[1]` = V3, `phi[2]` = V4, `gamma[1]` = V5, `gamma[2]` = V6) %>%
    gather(var, draw, -iter) %>%
    mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]',
                                        'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]'))) %>% 
    filter(iter > 1000) %>%
    ggplot() + geom_density(aes(draw)) + facet_wrap(~var, scales = 'free', ncol = 6, labeller = label_parsed) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 315, hjust = 0)) + 
    labs(x = NULL, y = NULL)
}


CHModel{
  
reps <- 80000
K <- 6
thin <- 10
burn <- 0.9

### TODO: Change hyper prior and initial values to reflect new model


draws <- list(list())
hyper <- list()
for(k in 1:K){
  hyper[[k]] <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(c(10, 10, 10, 10, 10, 10))), v = 6, scale = diag(1, 6))
  draws[[1]][[k]] <- list(mean = c(-5, -5, 0, 0, 0, 0), varInv = diag(10, 6))
}
draw[[1]]$pi <- rep(1/K, K)
hyper$alpha <- rep(1, K)
for(i in 1:N){
  draws[[i+1]] <- list(theta = c(-5, -5, 0, -0.1, 0.15, 0.05), k = sample(1:K, 1), pi = rep(1/K, K))
}
  
mixDraws <- mixtureMCMC(data, reps, draws, hyper, thin, K, 0.01)
  
saveRDS(mixDraws, 'mixN2000K6.RDS')
mixDraws <- readRDS(file = 'mixN2000K6.RDS')


### TODO: Clean up following analysis
  
  mapGroup <- NULL
  for(i in 1:N){
    kDraws <- mixDraws$draws[[i+1]]$k[(burn*reps/thin+1):(reps/thin)]
    mode <- which.max(table(c(1:K, kDraws)))
    mapGroup <- rbind(mapGroup, data.frame(group = c(mean(kDraws), mode), 
                                           method = c('mean', 'mode'),
                                           ID =  idSubset[i]))
  }
  ggplot(mapGroup) + geom_histogram(aes(group), binwidth = 0.1) + theme_bw() + facet_wrap(~method)
  
  muK <- NULL
  for(i in 1:K){
    mixDraws$draws[[1]][[i]]$mean[2:(reps/thin),] %>%
      cbind(iter = seq(2*thin, reps, thin)) %>%
      as.data.frame() %>%
      mutate(group = i)  -> temp
    colnames(temp) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2', 'iter', 'group')
    muK <- rbind(muK, temp)
  }
  
  muK %>%
    gather(var, draw, -iter, -group) %>%
    mutate(var = factor(var, levels = c('log_sigSq_eps', 'phi1', 'phi2', 'log_sigSq_eta', 'gamma1', 'gamma2'))) %>%
    filter(iter > 2000 & iter %% 20 == 0 ) %>%
    ggplot() + geom_line(aes(iter, draw)) + 
    facet_grid(var ~ group, scales = 'free') + theme_bw() + labs(title = 'mean') -> p1
  
  
  sdK <- NULL
  for(i in 1:K){
    mixDraws$draws[[1]][[i]]$varInv[,,2:(reps/thin)] %>%
      apply(3, function(x) sqrt(diag(solve(x)))) %>%
      t() %>%
      as.data.frame() %>%
      cbind(iter = seq(2*thin, reps, thin)) %>%
      mutate(group = i) -> temp
    colnames(temp) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2', 'iter', 'group')
    sdK <- rbind(sdK, temp)
  }
  
  sdK %>%
    gather(var, draw, -iter, -group) %>%
    mutate(var = factor(var, levels = c('log_sigSq_eps', 'phi1', 'phi2', 'log_sigSq_eta', 'gamma1', 'gamma2'))) %>%
    filter(iter > 2000 & iter %% 20 == 0 ) %>%
    ggplot() + geom_line(aes(iter, draw)) + 
    facet_grid(var ~ group, scales = 'free') + theme_bw() + labs(title = 'standard deviation') -> p2
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  
  for(i in 1:K){
    mixDraws$draws[[1]][[i]]$varInv[,,(reps/(2*thin)+1):(reps/thin)] %>% 
      apply(3, function(x) cov2cor(solve(x))) %>%
      t() %>%
      colMeans() %>%
      matrix(6) -> mat
    colnames(mat) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2')
    rownames(mat) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2') 
    corrplot::corrplot(mat) 
  }
  
  densities <- NULL
  support <- data.frame(seq(exp(-13), exp(-4.5), length.out = 1000),
                        seq(exp(-15), exp(-7.6), length.out = 1000),
                        seq(-0.5, 2, length.out = 1000),
                        seq(-1, 0.8, length.out = 1000),
                        seq(0, 2, length.out = 1000),
                        seq(-1, 0.5, length.out = 1000))
  vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'phi[1]', 'phi[2]', 'gamma[1]', 'gamma[2]') 
  index <- (burn*reps/thin+1):(reps/thin)
  for(k in 1:K){
    mat <- matrix(0, 1000, 6)
    for(i in seq_along(index)){
      meanvec <- mixDraws$draws[[1]][[k]]$mean[i,]
      sdvec <- sqrt(diag(solve(mixDraws$draws[[1]][[k]]$varInv[,,i])))
      w <- 0
      for(j in 1:N){
        w <- w + mixDraws$draws[[j+1]]$pi[i, k] / N
      }
      for(j in 1:2){
        mat[,j] <- mat[,j] + w * dlnorm(support[,j], meanvec[j], sdvec[j]) / length(index) 
      } 
      for(j in 3:6){
        mat[,j] <- mat[,j] + w * dnorm(support[,j], meanvec[j], sdvec[j]) / length(index) 
      }
    }
    densities <- rbind(densities,
                       data.frame(dens = c(mat),
                                  support = unlist(c(support)),
                                  var = rep(vars, rep(1000, 6)),
                                  group = k))
  }
  
  densities %>%
    mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]',
                                 'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]'))) -> densities
  densities %>% 
    group_by(var, support) %>%
    summarise(dens = sum(dens)) -> densMixMod
  
  densMixMod %>%
    ggplot(densMixMod) + geom_line(aes(support, dens)) +
    geom_line(data=densities, aes(support, dens, colour = factor(group))) +
    facet_wrap(~var, scales = 'free', ncol = 6, labeller = label_parsed) + 
    labs(x = NULL, y = NULL) + 
    theme_bw() + 
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))  -> p1
  
  
  ###TODO: Clean this for new model - estimation of prior
  
  means <- rep(0, 6*6)
  varinv <- matrix(0, 6*6, 6)
  linv <- matrix(0, 6*6, 6)
  pi <- rep(0, 6)
  for(i in 4001:8000){
    for(k in 1:6){
      means[(k-1)*6 + 1:6] <- means[(k-1)*6 + 1:6] + mixDraws$draws[[1]][[k]]$mean[i,] / 4000
      uinv <- chol(mixDraws$draws[[1]][[k]]$varInv[,,i])
      linv[(k-1)*6 + 1:6,] <- linv[(k-1)*6 + 1:6,] + uinv / 4000
      varinv[(k-1)*6 + 1:6,] <- varinv[(k-1)*6 + 1:6,] + t(uinv) %*% uinv / 4000
      for(n in 1:2000){
        pi[k] <- pi[k] + mixDraws$draws[[1+n]]$pi[i, k] / (2000 * 4000)
      }
    }
  }
  dets <- numeric(6)
  for(k in 1:6){
    dets[k] <- det(linv[(k-1)*6 + 1:6, ])  
  }
  priorMix <- list(mean = means, linv = linv, varInv = varinv, dets = dets, weights = pi)
  startingLam <- c(means, rep(c(diag(0.5, 6)), 6), rep(1, 6))
  
}

Independent{
  # Only used for comparision to E(theta | z) for the CH model
  lags <- 5
  reps <- 60000
  thin <- 10
  VAR <- FALSE
  
  error <- NULL
  for(lags in 1:5){
    for(i in 1:10){
      if(VAR){
        vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'rho', paste0('phi[', 1:(2*lags), ']'), paste0('gamma[', 1:(2*lags), ']'))
        draw <- c(-5, -5, rep(0, 1 + 4 * lags))
        hyper <- list(mean = c(-5, -5, rep(0, 1 + 4 * lags)), varInv = solve(diag(c(10, 10, 2.5, rep(10, 4 * lags)))))
        MCMC <- as.data.frame(indepMCMCVAR(data[[i]], reps, draw, hyper, lags = lags, thin = thin)$draws)
      } else {
        vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'rho', paste0('phi[', 1:lags, ']'), paste0('gamma[', 1:lags, ']'))
        draw <- c(-5, -5, rep(0, 1 + 2 * lags))
        hyper <- list(mean = c(-5, -5, rep(0, 1 + 2 * lags)), varInv = solve(diag(c(10, 10, 2.5, rep(10, 2 * lags)))))
        MCMC <- as.data.frame(indepMCMC(data[[i]], reps, draw, hyper, lags = lags, thin = thin)$draws)
      }
      
      colnames(MCMC) <- vars
      MCMC %>%
        mutate(iter = seq_along(rho) * thin) %>%
        filter(iter > reps/2) %>%
        mutate(`sigma[epsilon]^{2}` = exp(`sigma[epsilon]^{2}`),
               `sigma[eta]^{2}` = exp(`sigma[eta]^{2}`),
               rho = 2 / (1 + exp(-rho)) - 1) %>%
        colMeans() -> means
      
      T <- nrow(data[[i]])
      for(t in (T-100):(T-1)){
        a <- 0
        d <- 0
        for(j in 1:lags){
          if(VAR){
            a <- a  +  data[[i]][t - j, 1] * means[3 + j]  +  abs(data[[i]][t- j , 2]) * sign(data[[i]][t - j, 1]) * means[3 + lags + j]
            d <- d  +  data[[i]][t - j, 2] * means[3 + 2 * lags + j]  +  abs(data[[i]][t- j , 1]) * sign(data[[i]][t - j, 2]) * means[3 + 3 * lags + j]
          } else {
            a <- a  +  data[[i]][t - j, 1] * means[3 + j] 
            d <- d  +  data[[i]][t - j, 2] * means[3 + lags + j]
          }
        }
        error <- rbind(error, data.frame(a = (a - data[[i]][t, 1])^2,
                                         d = (d - data[[i]][t, 2])^2,
                                         lags = lags,
                                         var = VAR))
      }
      print(paste(lags, i, Sys.time()))  
    }
  }
  for(i in 1:10){
    T <- nrow(data[[i]])
    for(t in (T-100):(T-1)){
      error <- rbind(error, data.frame(a = data[[i]][t, 1]^2,
                                       d = data[[i]][t, 2]^2,
                                       lags = 0,
                                       var = VAR))
    }
  }
  error %>%
    gather(variable, SE, -lags, -var) %>%
    filter(is.finite(SE) & !is.na(SE)) %>%
    group_by(lags, variable, var) %>%
    summarise(MSE = mean(SE)) %>% 
    ggplot() + geom_line(aes(lags, MSE, colour = var)) + facet_wrap(~variable, scales = 'free') + ylim(0, 0.02)
  

}

# Combine results from the IH and CH model for this plot
  
posMeans %>%
  mutate(var = as.character(var),
        var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]',
                                      'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]'))) %>%
  filter((var == 'sigma[epsilon]^{2}' & mean > exp(-13) & mean < exp(-4.5)) |
           (var == 'sigma[eta]^{2}' & mean > exp(-15) & mean < exp(-7.6)) | 
           (var == 'phi[1]' & mean > -0.5 & mean < 2) |
           (var == 'phi[2]' & mean > -1 & mean < 0.8) | 
           (var == 'gamma[1]' & mean > 0 & mean < 2) | 
           (var == 'gamma[2]' & mean > -1 & mean < 0.5)) %>%
  ggplot() + geom_density(aes(mean)) + 
  geom_blank(data = densMixMod, aes(x = support)) + 
  facet_wrap(~var, scales = 'free', ncol = 6, labeller = label_parsed) +
  labs(x = NULL, y = NULL) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))-> p2
  
gridExtra::grid.arrange(p2, p1, ncol = 1)
  
  
}