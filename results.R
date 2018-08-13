library(tidyverse)
library(readr)

testid <- readRDS('idTest.RDS')

results <- list()
failed <- NULL
for(i in 1:1500){
  file <- paste0('results/car', testid[i], '.csv')
  file2 <- paste0('res2/car', testid[i], '.csv')
  if(file.exists(file) | file.exists(file2)){
 #   results[[i]] <- read_csv(file, col_types = cols())
  } else {
    failed <- c(failed, i)
  }
}

saveRDS(failed, 'failed.RDS')
failed <- readRDS('failed.RDS')

results <- bind_rows(results)
ggplot() + geom_bar(aes(floor((failed-1)/100)))

#write.csv(results, 'results.csv', row.names = FALSE)
#results <- readr::read_csv('results.csv')

saveRDS(results, 'results.RDS')
results <- readRDS('results.RDS')

results %>% 
  mutate(model = ifelse(model == 'Homogenous', 'Homogeneous', model)) -> results

# Average prediction Error

results[results$h %in% c(10, 20, 30) & results$inference %in% c('MCMC', NA),] %>%
  mutate(horizon = paste(h / 10, ifelse(h == 10, 'second', 'seconds'), 'ahead')) %>%
  group_by(S, horizon, model) %>%
  summarise(mean = mean(dist, na.rm = TRUE)) %>% 
  ggplot() + 
  geom_line(aes(S, mean, colour = model)) +
  facet_wrap(~horizon) + 
  labs(x = 'T', y = 'Mean Euclidean Error (metres)', colour = 'Model') + 
  theme_bw() +  
  theme(legend.position = 'bottom')


# Truncate to Time Series Only

results[results$h %in% c(10, 20, 30) & results$inference == 'MCMC',] %>%
  filter(!is.na(dist)) %>%
  mutate(horizon = paste(h / 10, ifelse(h == 10, 'second', 'seconds'), 'ahead')) %>%
  group_by(S, horizon, model) %>%
  summarise(mean = mean(dist, na.rm = TRUE)) %>% 
  ggplot() + 
  geom_line(aes(S, mean, colour = model)) +
  facet_wrap(~horizon) + 
  labs(x = 'T', y = 'Mean Euclidean Error (metres)', colour = 'Model') +
  theme_bw() +  
  theme(legend.position = 'bottom') 

# Rephrase the above two plots

results[results$inference %in% c('MCMC', NA), ] %>%
  group_by(h, model) %>%
  summarise(mean = mean(dist, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(h, mean, colour = model))

results[results$inference == 'MCMC', ] %>%
  filter(!is.na(logscore)) %>%
  mutate(S = ceiling((S-1)/100),
         S = paste0('observed up to ', S, '0 seconds')) %>%
  group_by(h, model, S) %>%
  summarise(mean = mean(dist, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(h, mean, colour = model)) + facet_wrap(~S)  + 
  labs(x = 'Forecast Horizon (100 milliseconds)', y = 'Mean Euclidean Error (metres)', colour = 'Model') +
  theme_bw() +  
  theme(legend.position = 'bottom') 


# Read in individual variances, plot against homogenous variances, and then split results by variance quintiles
sigma <- data.frame()
for(i in 1:15){
  temp <- read_csv(paste0('variances/group', i, '.csv'))
  sigma <- rbind(sigma, temp)
}
sigma$parameter <- as.factor(sigma$parameter)
# Get variance quantiles for each id, model, and parameter
sigma %>%
  group_by(parameter, model) %>%
  summarise(q1 = quantile(variance, 0.2),
            q2 = quantile(variance, 0.4),
            q3 = quantile(variance, 0.6), 
            q4 = quantile(variance, 0.8)) %>%
  right_join(sigma) %>%
  mutate(quantile = case_when(variance <= q1 ~ 1L,
                              variance <= q2 ~ 2L,
                              variance <= q3 ~ 3L,
                              variance <= q4 ~ 4L,
                              variance > q4 ~ 5L)) %>%
  select(-q1, -q2, -q3, -q4, -variance) %>%
  spread(parameter, quantile) -> varianceQuintiles

# Read in homogeneous variance means
homogDraws <- readRDS('homogMCMCN2000.RDS')$draws
homogMean = colMeans(exp(homogDraws[1001:5000, 1:2]))
homogMeanDf <- data.frame(parameter = factor(c('a', 'delta')),
                          posMean = colMeans(exp(homogDraws[1001:5000, 1:2])))

# Plot CH/IH variances against homogeneous
sigma %>%
  group_by(model, parameter) %>%
  filter(parameter == 'a' | variance < quantile(variance, 0.98)) %>%
  ggplot() + geom_density(aes(variance, colour = model)) + 
  geom_vline(data = homogMeanDf, aes(xintercept = posMean)) + 
  facet_wrap(~parameter, scales = 'free', ncol = 2) + 
  theme_bw() + 
  labs(x = 'Variance Posterior Mean', y = NULL, colour = 'Model')

# Percent of individual variances exceeded by homogeneous
sigma %>%
  left_join(homogMeanDf) %>%
  group_by(model, parameter) %>%
  summarise(percent = round(sum(variance < posMean) / n(), 3) * 100) %>%
  spread(model, percent)# %>%
  knitr::kable(format = 'latex')

# Compare each het. model relative to homog. for MCMC fits
results[results$h == 30 & results$inference == 'MCMC' & results$model %in% c('CH', 'IH', 'Homogeneous'), ] %>%
  select(-dist) %>%
  spread(model, logscore) %>%
  gather(model, logscore, CH, IH) %>%
  mutate(logscore = logscore - Homogeneous) %>%
  inner_join(varianceQuintiles) %>%
  group_by(a, delta, S, model) %>%
  summarise(logscore = median(logscore, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(S, logscore, colour = model)) + 
  geom_hline(aes(yintercept = 0), colour = 'black') + facet_grid(a ~ delta) + 
  labs(y = 'Median Difference in h = 30 Logscore (Heterogenous - Homogeneous)', x = 'T (100 milliseconds)', colour = 'Model') + 
  theme_bw() + 
  theme(legend.position = 'bottom')


# Plot updating logscores vs MCMC
results[!is.na(results$logscore) & 
          results$model %in% c('CH', 'IH') &
          results$h == 30,] %>%
  select(inference, logscore, S, id, model) %>%
  spread(inference, logscore) %>%
  mutate(Mixture = `UVB-Mixture` - MCMC,
         Single = `UVB-Single` - MCMC,
         S = ceiling((S - 100) / 50)) %>%
  gather(VB, logscore, Mixture, Single) %>%
  filter(S > 0) %>%
  ggplot() + geom_boxplot(aes(factor(S), logscore)) + 
  facet_grid(model~VB) + 
  labs(x = 'T (100 milliseconds)', y = 'Difference In Logscore (UVB - SVB)') + 
  theme_bw() +
  coord_cartesian(ylim = c(-5, 1.5)) + 
  scale_x_discrete(labels = c('110-150', '160-200', '210-250', '260-300', '310-350', '360-400', '410-450'))

# Compare CH UVB to Homogeneous MCMC
results[results$h == 30 & ((results$model %in% c('IH', 'CH') & results$inference %in% c('UVB-Single', 'UVB-Mixture'))
                            | results$model == 'Homogeneous'), ] %>%
  select(-dist) %>%
  mutate(type = paste(model, inference)) %>%
  select(-model, -inference) %>%
  spread(type, logscore) %>%
  gather(type, logscore, -h, -S, -id, -`Homogeneous MCMC`) %>%
  mutate(logscore = logscore - `Homogeneous MCMC`,
         model = substr(type, 1, 2)) %>%
  inner_join(varianceQuintiles) %>%
  group_by(a, delta, S, type) %>%
  filter(is.finite(logscore)) %>%
  summarise(logscore = median(logscore, na.rm = TRUE)) %>%
  filter(type %in% c('CH UVB-Mixture', 'IH UVB-Mixture')) %>%
  ggplot() + geom_line(aes(S, logscore, colour = type)) + 
  geom_hline(aes(yintercept = 0), colour = 'black') + facet_grid(a ~ delta) + 
  labs(y = 'Median Difference in h = 30 Logscore (Heterogeneous UVB - Homogeneous MCMC)', x = 'T (100 milliseconds)', colour = 'Model') + 
  theme_bw() + 
  theme(legend.position = 'bottom') 
