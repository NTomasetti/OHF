library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
source('MCMCFuns.R')

read.csv('carsAug.csv') %>%
  select(ID, relX, relY, relA, relV, relD) %>%
  ungroup() -> carsAug
colnames(carsAug) <- c('ID', 'x', 'a', 'y', 'v', 'delta')


# Select training vehicles from the vehicles that did not stop.
set.seed(1)
N <- 1000
carsAug %>%
  group_by(ID) %>%
  filter(min(relV) > 0) %>%
  ungroup() %>%
  .$ID %>%
  unique -> ids

idSubset <- sample(ids, N)
data <- list()
for(i in 1:N){
  carsAug %>%
    filter(ID == idSubset[i]) %>%
    mutate(d = relD - pi/2) %>% 
    select(relA , d, relV, relX, relY) %>%
    as.matrix() -> data[[i]]
}

saveRDS(data, 'MCMCData.RDS')
saveRDS(idSubset, 'trainingID.RDS')

# MCMC for the independent model
# This is used for model selection
Independent{
  
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
  for(i in 1:25){
    for(k in 1:2){
      for(lagsA in 1:4){
        for(lagsD in 1:4){
          T <- nrow(data[[i]])
          # Var = TRUE and FALSE have different parameters
          if(var[k]){
            vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'rho', paste0('phi[', 1:(2*lagsA), ']'), paste0('gamma[', 1:(2*lagsD), ']'))
            draw <- c(-5, -5, rep(0, 1 + 2 * (lagsA + lagsD)))
            hyper <- list(mean = c(-5, -5, rep(0, 1 + 2 * lagsA + 2 * lagsD)), varInv = solve(diag(c(10, 10, 2.5, rep(10, 2 * lagsA + 2 * lagsD)))))
          } else {
            vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'rho', paste0('phi[', 1:lagsA, ']'), paste0('gamma[', 1:lagsD, ']'))
            draw <- c(-5, -5, rep(0, 1 + lagsA + lagsD))
            hyper <- list(mean = c(-5, -5, rep(0, 1 + lagsA + lagsD)), varInv = solve(diag(c(10, 10, 2.5, rep(10, lagsA + lagsD)))))
          }
          MCMC <- as.data.frame(indepMCMC(data[[i]][1:(T-100),], reps, draw, hyper, lagsA = lagsA, lagsD = lagsD, thin = thin, VAR = var[k])$draws)
          
          # Transform the variance matrix from R^3 to sigma^2 and rho
          means <- c(0, 0, 0, colMeans(MCMC[1001:2000, 4:ncol(MCMC)]))

          for(t in (T-100):(T-10)){
            lags <- max(lagsA, lagsD)
            lagA <- data[[i]][-(0:lags) + t, 1]
            lagD <- data[[i]][-(0:lags) + t, 2] 
            dx <- 0
            dy <- 0
            v <- data[[i]][t, 3]
            for(h in 1:10){
              a <- 0
              d <- 0
              for(j in 1:lagsA){
                if(var[k]){
                  a <- a  +  lagA[j] * means[3 + j]  +  abs(lagD[j]) * sign(lagA[j]) * means[3 + lagsA + j]
                } else {
                  a <- a  +  lagA[j] * means[3 + j] 
                }
              }
              for(j in 1:lagsD){
                if(var[k]){
                  d <- d  +  lagD[j] * means[3 + 2 * lagsA + j]  +  abs(lagA[j]) * sign(lagD[j]) * means[3 + 2 * lagsA + lagsD + j]
                } else {
                  d <- d  +  lagD[j] * means[3 + lagsA + j]
                }
              }
              
              v <- v + a
              dx <- dx + v * cos(d + pi/2)
              dy <- dy + v * sin(d + pi/2)
              lagA <- c(a, lagA[1:(lags - 1)])
              lagD <- c(d, lagD[1:(lags - 1)])
            }
          error <- rbind(error, data.frame(error = sqrt((data[[i]][t+10, 4] - (data[[i]][t, 4] + dx))^2 + 
                                                          (data[[i]][t+10, 5] - (data[[i]][t, 5] + dy))^2),
                                           lagsA = lagsA,
                                           lagsD = lagsD,
                                           VAR = var[k],
                                           ID = i))
          }
        }
      }
    }
    print(paste0(i, ' ', Sys.time()))
  }
  
  # Compare MSE under each model
  error %>%
    group_by(lagsA, lagsD, VAR) %>%
    summarise(MSE = mean(error, na.rm = TRUE)) %>%
    arrange(MSE) %>%
    View()
  

  error %>%
    group_by(lagsA, lagsD, VAR) %>%
    summarise(MSE = mean(error, na.rm = TRUE)) %>%
    filter(VAR == TRUE) %>%
    ggplot() + geom_line(aes(lagsA, MSE, colour = factor(lagsD)))
  
  
  posMean %>%
    filter(VAR == TRUE & lagsA == 3 & lagsD == 2)  %>%
    mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'rho', paste0('phi[', 1:6, ']'),
                                        paste0('gamma[', 1:4, ']')))) %>%
    ggplot() + geom_density(aes(mean)) + facet_wrap(~var, scales = 'free', labeller = label_parsed)
  
  
}

# MCMC for the homogenous model
HomogenoousModel {
  reps <- 50000
  lagsA <- 4
  lagsD <- 4
  draws <- c(-5, -5, rep(0, 1 + lagsA + lagsD))
  hyper <- list(mean = c(-5, -5, rep(0, 1 + lagsA + lagsD)), varInv = solve(diag(c(10, 10, 2.5, rep(10, lagsA + lagsD)))))
  homogDraws <- homogMCMC(data, reps, draws, hyper, thin = 10, lagsA = lagsA, lagsD = lagsD, includeRho = FALSE, stepsize = 0.01)
  #saveRDS(homogDraws, 'homogMCMCN3000.RDS')
  
  homogDraws$draws %>%
    as.data.frame() %>%
    mutate(iter = seq_along(V1) * thin,
           V1 = exp(V1), 
           V2 = exp(V2),
           V3 = 2 / (1 + exp(-V3)) - 1) %>%
    rename(`sigma[epsilon]^{2}` = V1, `sigma[eta]^{2}` = V2, `rho` = V3, `phi[1]` = V4, `phi[2]` = V5, `phi[3]` = V6, `gamma[1]` = V7, `gamma[2]` = V8) %>%
    gather(var, draw, -iter) %>%
    mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]', 'phi[3]',
                                        'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]', 'rho'))) %>% 
    filter(iter > reps / 2) %>%
    ggplot() + geom_line(aes(iter, draw)) + 
    facet_wrap(~var, scales = 'free', labeller = label_parsed) +
    theme_bw()
    
  homogDraws$draws %>%
    as.data.frame() %>%
    mutate(iter = seq_along(V1) * thin,
           V1 = exp(V1), 
           V2 = exp(V2),
           V3 = 2 / (1 + exp(-V3)) - 1) %>%
    rename(`sigma[epsilon]^{2}` = V1, `sigma[eta]^{2}` = V2, `rho` = V3, `phi[1]` = V4, `phi[2]` = V5, `phi[3]` = V6, `gamma[1]` = V7, `gamma[2]` = V8) %>%
    gather(var, draw, -iter) %>%
    mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]', 'phi[3]',
                                        'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]', 'rho'))) %>% 
    filter(iter > reps / 2) %>%
    ggplot() + geom_density(aes(draw)) + facet_wrap(~var, scales = 'free', labeller = label_parsed) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 315, hjust = 0)) + 
    labs(x = NULL, y = NULL)
}

# MCMC for the clustered heterogeneous model
CHModel{
  
reps <- 50000
K <- 3
thin <- 10
burn <- 0.9
lagsA <- 3
lagsD <- 2

draws <- list(list())
hyper <- list()
for(k in 1:K){
  hyper[[k]] <- list(mean = c(-5, -5, rep(0, 1 + lagsA + lagsD)),
                     varInv = solve(diag(c(10, 10, 2.5, rep(10, lagsA + lagsD)))),
                     v = 3 + lagsA + lagsD, 
                     scale = diag(1, 3 + lagsA + lagsD))
  draws[[1]][[k]] <- list(mean = c(-5, -5, rep(0, 1 + lagsA + lagsD)), varInv = diag(10, 3 + lagsA + lagsD))
}
hyper$alpha <- rep(1, K)
for(i in 1:N){
  draws[[i+1]] <- list(theta = c(-5, -5, rep(0, 1 + lagsA + lagsD)), k = sample(1:K, 1), pi = rep(1/K, K))
}
  
mixDraws <- mixtureMCMC(data, reps, draw, hyper, thin , K, lagsA, lagsD)
  
saveRDS(mixDraws, paste0('mixN', N, 'K', K, '.RDS'))
mixDraws <- readRDS(file = paste0('mixN', N, 'K', K, '.RDS'))


# Weights
mixDraws$draws[[1]]$pi %>%
  as.data.frame() %>%
  mutate(iter = seq_along(V1) * thin) %>%
  filter(iter > burn*reps) %>%
  gather(pi, draw, -iter) %>%
  ggplot() + geom_line(aes(iter, draw)) + 
  facet_wrap(~pi, scales = 'free')

mixDraws$draws[[1]]$pi %>%
  as.data.frame() %>%
  mutate(iter = seq_along(V1) * thin) %>%
  filter(iter >  burn*reps) %>%
  gather(pi, draw, -iter) %>%
  ggplot() + geom_density(aes(draw)) + 
  facet_wrap(~pi, scales = 'free')
  
muK <- NULL
for(i in 1:K){
  mixDraws$draws[[1]][[i]]$mean %>%
    as.data.frame() %>%
    rename(`sigma[epsilon]^{2}` = V1, `sigma[eta]^{2}` = V2, `rho` = V3, `phi[1]` = V4, `phi[2]` = V5, `phi[3]` = V6, `gamma[1]` = V7, `gamma[2]` = V8) %>%
    mutate(iter = seq_along(rho) * thin,
           group = i) %>%
    filter(iter >  burn*reps) %>%
    rbind(muK) -> muK
}
  
muK %>%
  gather(var, draw, -iter, -group) %>%
  mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]', 'phi[3]',
                                      'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]', 'rho'))) %>% 
  ggplot() + geom_line(aes(iter, draw)) + 
  facet_grid(var ~ group, scales = 'free') + theme_bw() + labs(title = 'mean') -> p1
  
  
sdK <- NULL
for(i in 1:K){
  mixDraws$draws[[1]][[i]]$varInv[,,2:(reps/thin)] %>%
      apply(3, function(x) sqrt(diag(solve(x)))) %>%
      t() %>%
      as.data.frame() %>%
    rename(`sigma[epsilon]^{2}` = V1, `sigma[eta]^{2}` = V2, `rho` = V3, `phi[1]` = V4, `phi[2]` = V5, `phi[3]` = V6, `gamma[1]` = V7, `gamma[2]` = V8) %>%
    mutate(iter = seq_along(rho) * thin,
           group = i) %>%
    filter(iter > burn * reps) %>%
    rbind(sdK) -> sdK
}
  
sdK %>%
    gather(var, draw, -iter, -group) %>%
    mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]', 'phi[3]',
                                        'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]', 'rho'))) %>% 
    ggplot() + geom_line(aes(iter, draw)) + 
    facet_grid(var ~ group, scales = 'free') + theme_bw() + labs(title = 'standard deviation') -> p2
  
gridExtra::grid.arrange(p1, p2, ncol = 2)
  
densities <- NULL
support <- data.frame(seq(exp(-13), exp(-6.5), length.out = 1000),
                      seq(exp(-15), exp(-7.5), length.out = 1000),
                      seq(-1, 1, length.out = 1000),
                      seq(-2, 2, length.out = 1000),
                      seq(-2, 2, length.out = 1000),
                      seq(-2, 2, length.out = 1000),
                      seq(-2, 2, length.out = 1000),
                      seq(-2, 2, length.out = 1000))
vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'rho', 'phi[1]', 'phi[2]', 'phi[3]', 'gamma[1]', 'gamma[2]') 
index <- (reps/thin+1):(reps/thin)

rhoDens <- function(x, mu, sd){
  abs(1 / (1 + x) + 1 / ( 1- x)) * dnorm(log((x+1)/(1-x)), mu, sd)
}

for(k in 1:K){
  mat <- matrix(0, 1000, 8)
  for(i in seq_along(index)){
    meanvec <- mixDraws$draws[[1]][[k]]$mean[index[i],]
    sdvec <- sqrt(diag(solve(mixDraws$draws[[1]][[k]]$varInv[,,index[i]])))
    w <- mixDraws$draws[[1]]$pi[index[i], k]
    for(j in 1:2){
      mat[,j] <- mat[,j] + w * dlnorm(support[,j], meanvec[j], sdvec[j]) / length(index) 
    } 
    mat[, 3] <- mat[, 3] + w * rhoDens(support[,3], meanvec[3], sdvec[3]) / length(index)
    for(j in 4:8){
      mat[,j] <- mat[,j] + w * dnorm(support[,j], meanvec[j], sdvec[j]) / length(index) 
    }
  }
  densities <- rbind(densities,
                     data.frame(dens = c(mat),
                                support = unlist(c(support)),
                                var = rep(vars, rep(1000, 8)),
                                group = k))
}
densities %>%
  mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]', 'phi[3]',
                                 'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]', 'rho'))) -> densities
densities %>% 
  group_by(var, support) %>%
  summarise(dens = sum(dens)) -> densMixMod
  
ggplot(densMixMod) + geom_line(aes(support, dens)) +
  geom_line(data=densities, aes(support, dens, colour = factor(group))) +
  facet_wrap(~var, scales = 'free', labeller = label_parsed) + 
  labs(x = NULL, y = NULL) + 
  theme_bw() + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))

means <- rep(0, K*8)
varinv <- matrix(0, K*8, 8)
linv <- matrix(0, K*8, 8)
pi <- rep(0, K)
for(i in seq_along(index)){
  for(k in 1:K){
    means[(k-1)*8 + 1:8] <- means[(k-1)*8 + 1:8] + mixDraws$draws[[1]][[k]]$mean[index[i],] / length(index)
    uinv <- chol(mixDraws$draws[[1]][[k]]$varInv[,,index[i]])
    linv[(k-1)*8 + 1:8,] <- linv[(k-1)*8 + 1:8,] + t(uinv) / length(index)
    varinv[(k-1)*8 + 1:8,] <- varinv[(k-1)*8 + 1:8,] + uinv %*% t(uinv) / length(index)
    pi[k] <- pi[k] + mixDraws$draws[[1]]$pi[index[i],k] / length(index)
  }
}
dets <- numeric(K)
for(k in 1:K){
  dets[k] <- det(linv[(k-1)*8 + 1:8, ])  
}
priorMix <- list(mean = means, linv = linv, varInv = varinv, dets = dets, weights = pi)
startingLam <- c(means, rep(c(diag(0.5, 8)), K), rep(1, K))
  
}


