library(tidyverse)
library(Rcpp)
library(RcppArmadillo)

# Load data.
cars1 = read.table('vehicle-trajectory-data/trajectories1.txt')

cars2 = read.table('vehicle-trajectory-data/trajectories2.txt')
cars2[,1] = cars2[,1] + 10000 # IDs restart at one for each segment, so add 10000/20000 so IDs are unique for the combined set
cars2[,15] = cars2[,15] + 10000

cars3 = read.table('vehicle-trajectory-data/trajectories3.txt')
cars3[,1] = cars3[,1] + 20000
cars3[,15] = cars3[,15] + 20000

cars = rbind(cars1, cars2, cars3)
rm(cars1, cars2, cars3)

colnames(cars) = c('ID', 'frame', 'totalFrames', 'time', 'x', 'y', 
                   'globalX', 'globalY', 'length', 'width', 'class',
                   'veloc', 'accel', 'lane', 'proceeding', 'following', 
                   'spacing', 'headway')

# Example data picture
cars %>%
  filter(ID < 100) %>%
  ggplot() + geom_path(aes(x, y, group = ID)) +
  labs(y = NULL, x = 'Lane') + 
  theme_bw() + 
  scale_x_continuous(labels = c('Lane 1', 'Lane 2', 'Lane 3', 'Lane 4', 'Lane 5', 'Entry Lane'),
                     breaks = c(6, 18, 30, 40, 50, 65))

# Initial operations on data lanes, and calculating distance
cars %>%
  group_by(ID) %>%
  summarise(medLane = median(lane),
            changed = any(lane != medLane),
            enterExit = any(lane > 5),
            startLane = head(lane, 1)) %>%
  ungroup() %>%
  right_join(cars, by = "ID") -> cars

cars %>%
  filter(enterExit == FALSE) %>%
  group_by(ID) %>%
  mutate(time = time - min(time),
         n = seq_along(frame), 
         xlag = ifelse(n == 1, 0, lag(x)), 
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x-xlag)^2 + (y-ylag)^2),
         delta = atan2(y - ylag, x - xlag),
         dist = cumsum(v)) -> cars

# Define spline degree and extract 100 vehicles per lane that did not change lane to estimate spline models
degree <- 50

cars %>%
  filter(changed == FALSE) %>%
  group_by(lane) %>%
  filter(ID %in% head(unique(ID), 100)) %>%
  .$ID %>%
  unique() -> splinesID

# Input pos = (x, y) and centre, a dataframe of distance / spline xhat / spline yhat to transform coordinates
# The function will find the element of the xhat / yhat combo of centre which minimises the Euc. distance to the input (x, y).
transformCoord <- function(pos, centre){
  x <- pos[1]
  y <- pos[2]
  centre <- centre[abs(centre$yhat - y) < 15,] 
  centre$dist <- sqrt((x - centre$xhat)^2 + (y - centre$yhat)^2)
  closest <- centre[which.min(centre$dist),]
  relX <- sign(x - closest$xhat) * closest$dist
  relY <- closest$d
  return(c(relX, relY))
}

# Loop through each lane to:
# 1) Fit spline models
# 2) Evaluate spline on a fine grid of distances
# 3) Loop through each vehicle in the lane to calculate relative coordinates
# Result: Relcoord dataframe with relX / relY / time / ID
relCoord <- NULL
skip <- NULL
for(i in 1:5){
  cars %>%
    filter(lane == i  & ID %in% splinesID) -> carSubset
  modX <- smooth.spline(carSubset$dist, carSubset$x, df = degree)
  modY <- smooth.spline(carSubset$dist, carSubset$y, df = degree)

  cars %>%
    filter(startLane == i & !ID %in%splinesID) %>% 
    ungroup() -> toPredict
  
  centrePath <- data.frame(d = seq(min(toPredict$dist), max(toPredict$dist), length.out = 1000000))
  centrePath$xhat <- predict(modX, centrePath$d)$y
  centrePath$yhat <- predict(modY, centrePath$d)$y
  
  for(j in seq_along(unique(toPredict$ID))){
    carJ <- filter(toPredict, ID == unique(toPredict$ID)[j])
    rel <- select(carJ, x, y) %>%
      apply(1, function(row) transformCoord(c(row[1], row[2]), centrePath)) %>% 
      t() %>% as.data.frame()
    if(nrow(rel) == 1){
      skip <- c(skip, carJ$ID[1])
      next
    }
    
    colnames(rel) <- c('relX', 'relY')
    rel$ID <- carJ$ID
    rel$time <- carJ$time
    
    relCoord <- rbind(relCoord, rel)
    print(paste(i, j))
  }
}
mPerFoot <- 0.3048

# Combine relcoord with the original cars dataset
# Calculate relative velocity / angle / acceleration / change in angle

cars %>% 
  filter(!ID %in% c(skip, splinesID)) %>%
  left_join(relCoord) %>%
  group_by(ID) %>%
  mutate(n = seq_along(time),
         relX = relX * mPerFoot,
         relY = relY * mPerFoot,
         x = x * mPerFoot,
         y = y * mPerFoot,
         relV = sqrt((relX - lag(relX))^2 + (relY - lag(relY))^2),
         relA = relV - lag(relV),
         relD = atan2(relY - lag(relY), relX - lag(relX)),
         changeD = relD - lag(relD)) %>%
  filter(n > 2) %>%
  ungroup() %>% 
  select(ID, relX, relY, relV, relA, relD, changeD, time, changed, lane, startLane, x, y, time) -> carsAug

# A c++ funciton to loop through delta and find cases where relV == 0, and apply d_t = d_{t-1}

cppFunction(depends = "RcppArmadillo",
            'arma::vec lagDelta(arma::mat data) {
              int T = data.n_rows;
              arma::vec out(T);
              out(0) = data(0, 5);
              for(int i = 1; i < T; ++i){
                if(data(i, 0) == data(i-1, 0) & data(i, 3) == 0){
                  out(i) = out(i-1);
                } else {
                  out(i) = data(i, 5);
                }
              }
              return out;
            }')

carsAug$relD = lagDelta(as.matrix(carsAug))


write.csv(carsAug, 'carsAug.csv', row.names = FALSE)

# This will plot the time series of v / d / a / change in d and associated pacfs for a sample of vehicles.
# Even with the d_t = d_{t-1} transform for vehicles that have stopped, stationary (y_t = y_{t-1}) vehicles can have tiny movements in x that lead to huge spikes in delta to 0 or pi.
# This may not matter too much, if V is small the angle does not really matter.
# Drop the filter / group / ungroup lines to include vehicles that stopped.

carsAug %>%
  group_by(ID) %>% 
  filter(min(relV) > 0) %>%
  ungroup() %>%
  .$ID %>%
  unique() %>%
  sample(5) -> plotID

carsAug %>%
  filter(ID %in% plotID) -> pacfData

pacfData %>%
  select(ID, time, relD, relA, relV) %>%
  group_by(ID) %>%
  mutate(time = time - min(time)) %>%
  rename(a = relA, v = relV, del = relD) %>%
  ungroup() %>%
  mutate(ID = case_when(ID == plotID[1] ~ 'Vehicle 1', 
                        ID == plotID[2] ~ 'Vehicle 2',
                        ID == plotID[3] ~ 'Vehicle 3',
                        ID == plotID[4] ~ 'Vehicle 4',
                        ID == plotID[5] ~ 'Vehicle 5'),
         ID = factor(ID)) %>%
  ungroup() %>%
  gather(var, value, -ID, -time) %>%
  ggplot() + geom_line(aes(time, value)) + 
  facet_grid(var ~ ID, scales = 'free') + 
  labs(x = 'T', y = 'Value') + theme_bw()


pacfs <- data.frame()
for(i in 1:5){
  id <- unique(pacfData$ID)[i]
  
  pacfData %>%
    filter(ID == id) %>%
    .$relD %>%
    pacf(plot = FALSE) %>%
    with(data.frame(lag, acf)) %>%
    mutate(ID = paste0('Vehicle ', i),
           var = 'del') -> tmp
  pacfs <- rbind(pacfs, tmp)
  
  pacfData %>%
    filter(ID == id) %>%
    .$relA %>%
    pacf(plot = FALSE) %>%
    with(data.frame(lag, acf)) %>%
    mutate(ID = paste0('Vehicle ', i),
           var = 'a') -> tmp
  pacfs <- rbind(pacfs, tmp)
  
  pacfData %>%
    filter(ID == id) %>%
    .$relV %>%
    pacf(plot = FALSE) %>%
    with(data.frame(lag, acf)) %>%
    mutate(ID = paste0('Vehicle ', i),
           var = 'v') -> tmp
  pacfs <- rbind(pacfs, tmp)
  
  pacfData %>%
    filter(ID == id) %>%
    .$changeD %>%
    pacf(plot = FALSE) %>%
    with(data.frame(lag, acf)) %>%
    mutate(ID = paste0('Vehicle ', i),
           var = 'changeD') -> tmp
  pacfs <- rbind(pacfs, tmp)
  
}

ggplot(pacfs) + geom_bar(aes(lag, acf), stat = 'identity') + 
  facet_grid(var~ID) + labs(x = 'Lag', y = 'Partial Autocorrelation') + 
  theme_bw()

# Calculate average distance between lanes for lane change prediction

d <- seq(min(cars$dist), max(cars$dist), length.out = 1000000)
centrePath <- matrix(d, ncol = 1)
for(i in 1:5){
  cars %>%
    filter(lane == i  & ID %in% splinesID) -> carSubset
  modX <- smooth.spline(carSubset$dist, carSubset$x, df = degree)
  modY <- smooth.spline(carSubset$dist, carSubset$y, df = degree)
  centrePath <- cbind(centrePath, predict(modX, d)$y)
  centrePath <- cbind(centrePath, predict(modY, d)$y)
}
centrePath <- as.data.frame(centrePath)
colnames(centrePath) <- c('d', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'x5', 'y5')

centrePath <- mutate(centrePath,
                     g1 = sqrt((x1 - x2)^2 + (y1 - y2)^2),
                     g2 = sqrt((x2 - x3)^2 + (y2 - y3)^2),
                     g3 = sqrt((x3 - x4)^2 + (y3 - y4)^2),
                     g4 = sqrt((x4 - x5)^2 + (y4 - y5)^2))
laneSize <- mean(unlist(select(centrePath, g1, g2, g3, g4))) * mPerFoot

