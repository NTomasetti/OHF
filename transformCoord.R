library(tidyverse)

# Load data
cars1 = read.table('vehicle-trajectory-data/trajectories1.txt')

cars2 = read.table('vehicle-trajectory-data/trajectories2.txt')
cars2[,1] = cars2[,1] + 10000
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

cars %>%
  filter(ID %in% head(unique(.$ID), 10)) %>%
  ggplot() + geom_path(aes(x, y, group = ID)) +
  labs(y = NULL, x = 'Lane') + 
  theme_bw() + 
  scale_x_continuous(labels = c('Lane 1', 'Lane 2', 'Lane 3', 'Lane 4', 'Lane 5', 'Entry Lane'),
                     breaks = c(6, 18, 30, 40, 50, 65))

#Operations on data
cars %>%
  mutate(time = time - min(time)) -> cars
cars %>%
  group_by(ID) %>%
  summarise(medLane = median(lane),
            changed = any(lane != medLane),
            enterExit = any(lane > 5)) %>%
  ungroup() %>%
  right_join(cars, by = "ID") -> cars

cars %>% 
  group_by(ID) %>%
  mutate(n = seq_along(time)) %>%
  filter(n == 1) %>%
  select(ID, lane) %>%
  rename(startLane = lane) %>%
  left_join(cars) -> cars

cars %>%
  group_by(ID) %>%
  filter(enterExit == FALSE) %>%
  mutate(n = seq_along(frame), 
         xlag = ifelse(n == 1, 0, lag(x)), 
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x-xlag)^2 + (y-ylag)^2),
         delta = atan2(y - ylag, x - xlag),
         dist = cumsum(v)) -> cars



degree <- 50

cars %>%
  filter(changed == FALSE) %>%
  group_by(lane) %>%
  filter(ID %in% head(unique(ID), 100)) %>%
  .$ID %>%
  unique() -> splinesID

transformCoord <- function(pos, centre){
  x <- pos[1]
  y <- pos[2]
  centre[abs(centre$yhat - y) < 15,] -> centre
  centre$dist <- sqrt((x - centre$xhat)^2 + (y - centre$yhat)^2)
  closest <- centre[which.min(centre$dist),]
  relX <- sign(x - closest$xhat) * closest$dist
  relY <- closest$d
  
  return(c(relX, relY))
}

relCoord <- NULL
for(i in 1:5){
  cars %>%
    filter(lane == i  & ID %in% splinesID) -> carSubset
  modX <- smooth.spline(carSubset$dist, carSubset$x, df = degree)
  modY <- smooth.spline(carSubset$dist, carSubset$y, df = degree)

  cars %>%
    filter(startLane == i & !ID %in%splinesID) %>% 
    ungroup() -> toPredict
  
  centrePath <- data.frame(d = seq(min(toPredict$dist), max(toPredict$dist), length.out = 25000))
  centrePath$xhat <- predict(modX, centrePath$d)$y
  centrePath$yhat <- predict(modY, centrePath$d)$y
  
  for(j in seq_along(unique(toPredict$ID))){
    carJ <- filter(toPredict, ID == unique(toPredict$ID)[j])
    rel <- select(carJ, x, y) %>%
      apply(1, function(row) transformCoord(c(row[1], row[2]), centrePath)) %>% 
      t() %>% as.data.frame()
    colnames(rel) <- c('relX', 'relY')
    rel$ID <- unique(toPredict$ID)[j]
    rel$time <- carJ$time
    
    relCoord <- rbind(relCoord, rel)
    print(paste(i, j))
  }
}

cars %>% 
  filter(!ID %in%splinesID) %>%
  left_join(relCoord) %>%
  group_by(ID) %>%
  mutate(n = seq_along(time),
         lagRelX = ifelse(n == 1, 0, lag(relX)),
         lagRelY = ifelse(n == 1, 0, lag(relY)),
         relV = sqrt((relX - lagRelX)^2 + (relY - lagRelY)^2),
         relA = relV - ifelse(n <= 2, 0, lag(relV)),
         relD = atan2(relY - lagRelY, relX - lagRelX)) %>%
  filter(n > 2) %>%
  ungroup() %>% 
  select(ID, changed, lane, startLane, relX, relY, relV, relA, relD, x, y, time) -> carsAug

write.csv(carsAug, 'carsAug.csv', row.names = FALSE)