library(tidyverse)
library(viridis)



### Preparation code #########################################################

# read the file with simulation results you are interested in plotting
AICall <- read.csv(file="results_simulations/simulationCluster.csv", sep =" ")

# specify the title of plots (with respect to the simulation type you've chosen)
titleOfPlots <- ""
#titleOfPlots <- "Type F"



### FIG 1 ####################################################################
library(plot.matrix)

teoret<-matrix(0, nrow=50, ncol=50)

for(i in 1:50){
  for(j in 1:50){
    teoret[i,j]<-ifelse(i<j,(i/j)^0.5, NA)
  }}

par(mar=c(5,5,5,5))
plot(teoret, main="Expected (theoretical) correlation between spatial lags 
     for different knn (used in W)", xlab="knn in W1", ylab="knn in W2")



### FIG 2 ####################################################################

# Flowchart prepared in Canva software



### FIG 3 ####################################################################

# simulation design - small sample, middle size sample, big sample, colored by roa
# code stored in the "Structure_Spatial_Plot with general sample structure in three sizes.R" file



### FIG 4 ####################################################################

# Nominal (raw) values of AIC due to incremental increase of knn in W
# y = AIC, x = knn, small medium big sample

head(AICall)
AICall <- AICall[-1,] #remove first row

sampleLabs <- c("200 points", "500 points", "1000 points")
names(sampleLabs) <- c(200, 500, 1000)

summary(AICall)

ggplot(data = AICall, aes(x = knn, y=AIC, group = seed)) + geom_line(aes(color = seed))+
  facet_wrap(~obsSample, scales = "free",  labeller = labeller(obsSample = sampleLabs))+
  scale_color_viridis_c()+
  labs(colour="seed", title = titleOfPlots)+
  theme_minimal()+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))

#export 1250, 300


### FIG 5 ####################################################################

# Relative values of AIC against relative values of knn
# AIC standarized (to 0-1 scale) on y, on x - rank, where AIC reaches minimum =0, and then consecutive +1
# nor egression lines. U-shaped relation

AICall$run <- NA
AICall$run <- rep(1:60, each = 49)
summary(AICall)

data01 <- AICall %>% group_by(obsSample, seed) %>% summarise(AICmin = min(AIC), AICmax = max(AIC))
data01

#standarizing AIC to 0-1 interval 
# (AIC-AICmin)/(AICmax-AICmin)
data02 <- merge(AICall, data01) %>% mutate(AICrel = (AIC - AICmin)/(AICmax -AICmin)) %>% arrange(run)
data02

dataxx <- cbind(data02$knn, data02$AIC, data02$AICrel, data02$rho, data02$obsSample, data02$run, data02$seed)
colnames(dataxx) <- c('knn', 'AIC',  'AICrel', 'rho', 'n', 'run', 'seed')
dataxx <- as.data.frame(dataxx)
summary(dataxx)


dataRank <- dataxx %>% group_by(run) %>% mutate(rank = row_number() - which.min(AIC))

ggplot(data = dataRank, aes(x = rank, y=AICrel, group = seed)) + geom_point(aes(color = seed))+
  facet_wrap(~n, scales = "free",  labeller = labeller(n = sampleLabs))+
  scale_color_viridis_c()+
  labs(colour="seed", x = 'knn.rel', y = 'AICst', title = titleOfPlots)+
  theme_minimal()+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))

#export 1250, 300
### export 950, 400


### FIG 6 ####################################################################

# Values of relative rho (rho/rho(minAIC)) against relative knn (knn=0 for minAIC)

data01 <- AICall %>% group_by(obsSample, seed) %>% summarise(AICmin = min(AIC), AICmax = max(AIC))
data01

#standarizing AIC to 0-1 interval 
# (AIC-AICmin)/(AICmax-AICmin)
data02 <- merge(AICall, data01) %>% mutate(AICrel = (AIC - AICmin)/(AICmax -AICmin)) %>% arrange(run)
data02


dataxx <- cbind(data02$knn, data02$AIC, data02$AICrel, data02$rho, data02$obsSample, data02$run, data02$seed)
colnames(dataxx) <- c('knn', 'AIC',  'AICrel', 'rho', 'n', 'run', 'seed')
dataxx <- as.data.frame(dataxx)
summary(dataxx)


dataRank <- dataxx %>% group_by(run) %>% mutate(rank = row_number() - which.min(AIC))


minRhoByRank <- dataRank %>% group_by(run) %>% mutate(rho.minAIC = ifelse(rank==0, rho, 0)) %>% 
  summarise(run = mean(run), rho.min.AIC = sum(rho.minAIC))

dataRank2 <- merge(dataRank, minRhoByRank)  %>% mutate(rho.rel = rho/rho.min.AIC) %>% arrange(run)

#Figure 5 rho rel vs rank free scale

ggplot(data = dataRank2, aes(x = rank, y=rho.rel, group = seed)) + geom_point(aes(color = seed))+
  facet_wrap(~n, scales = "free",  labeller = labeller(n = sampleLabs))+
  scale_color_viridis_c()+
  labs(colour="seed", y = 'rho/rho(minAIC)', x = 'knn.rel', title = titleOfPlots)+
  theme_minimal()+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))

#export 950, 400



### FIG 7 ####################################################################

# Values of relative total impact (total/total(minAIC)) against relative knn (knn=0 for minAIC)

data01 <- AICall %>% group_by(obsSample, seed) %>% summarise(AICmin = min(AIC), AICmax = max(AIC))
data01

#standarizing AIC to 0-1 interval 
# (AIC-AICmin)/(AICmax-AICmin)
data02 <- merge(AICall, data01) %>% mutate(AICrel = (AIC - AICmin)/(AICmax -AICmin)) %>% arrange(run)
data02


dataxx <- cbind(data02$knn, data02$AIC, data02$AICrel, data02$total, data02$obsSample, data02$run, data02$seed)
colnames(dataxx) <- c('knn', 'AIC',  'AICrel', 'total', 'n', 'run', 'seed')
dataxx <- as.data.frame(dataxx)
summary(dataxx)


dataRank <- dataxx %>% group_by(run) %>% mutate(rank = row_number() - which.min(AIC))


minTotalByRank <- dataRank %>% group_by(run) %>% mutate(total.minAIC = ifelse(rank==0, total, 0)) %>% 
  summarise(run = mean(run), total.min.AIC = sum(total.minAIC))

dataRank2 <- merge(dataRank, minTotalByRank)  %>% mutate(total.rel = total/total.min.AIC) %>% arrange(run)


ggplot(data = dataRank2, aes(x = rank, y=total.rel, group = seed)) + geom_point(aes(color = seed))+
  facet_wrap(~n, scales = "free",  labeller = labeller(n = sampleLabs))+
  scale_color_viridis_c()+
  labs(colour="seed", y = 'total/total(minAIC)', x = 'knn.rel', title = paste0(titleOfPlots, " total impact (dist)"))+
  theme_minimal()+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))

#export 950, 400


### FIG 8 ####################################################################

# Moran for dependent variable with regard to incremental increasing knn

sampleLabs <- c("200 points", "500 points", "1000 points")
names(sampleLabs) <- c(200, 500, 1000)

summary(AICall)

ggplot(data = AICall, aes(x = knn, y=moran, group = seed)) + geom_line(aes(color = seed))+
  facet_wrap(~obsSample, scales = "free",  labeller = labeller(obsSample = sampleLabs))+
  scale_color_viridis_c()+
  labs(colour="seed")+
  theme_minimal()+
  labs(y = "Moran's I statistic", title = titleOfPlots)+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))

#export 950, 400



### FIG 9 ####################################################################

# Variogram-like statistic against knn 

#semiVariogram analysis for type A simulation
dane<-read.csv("results_simulations/semiVarianceClusters.csv", sep=" ")
dane
dane<-dane[-1,] #remove empty first row
dane$seed <- seq(1,1000,50)
dane$obs <- c(rep(200,20), rep(500,20), rep(1000,20))

dane<-as.data.frame(dane)
dane$knn <- c(2:50)
dane<-dane[,-1] #remove empty first column

head(dane)

dane <- dane %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "knn",
    names_prefix = "V",
    values_to = "semiVar",
    values_drop_na = TRUE
  ) 

dane <- dane %>% 
  mutate(knn = as.numeric(knn)) %>% 
  mutate(seed = as.numeric(seed))

ggplot(data = dane, aes(x = knn, y=semiVar, group = seed)) + geom_line(aes(color = seed))+
  facet_wrap(~obs, scales = "free",  labeller = labeller(obs = sampleLabs))+
  scale_color_viridis_c()+
  labs(colour="seed")+
  theme_minimal()+
  labs(y = "Variogram-like statistic")+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))

#export 950, 400


### FIG 10 ####################################################################

# Empirical data sample, size 200, 500, 1000 and 10k
# AIC vs knn

AICempi <- read.csv(file="results_simulations/resultsEmpiricalSimulation.csv", sep =";", dec = ",")
head(AICempi)

sampleLabs2 <- c("200 points", "500 points", "1000 points", "10000 points" )
names(sampleLabs2) <- c(200, 500, 1000, 10000)

library(viridis)

ggplot(data = AICempi, aes(x = knn, y=AIC, group = minInterval)) + geom_line(aes(color = minInterval))+
  facet_wrap(~obsSample, scales = "free",  labeller = labeller(obsSample = sampleLabs2))+
  scale_color_viridis_c()+
  labs(colour="subset ID")+
  theme_minimal()+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))
