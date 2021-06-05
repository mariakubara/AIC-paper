library(tidyverse)

# Generating datasets

# setting the seed for reproducibility 
set.seed(601) 
#set.seed(12345678)

datax <- matrix(NA, ncol =16)
colnames(datax) <- c("xi", "yi", "powLub", "roa_geo", "subsample", "sector", "roa_sec","dist", "roa_dist" ,
                     "roa_param", "roa", "agri", "prod", "constr", "serv", "sample" )

numObs<-list(c(6,25,2,25),c(15,25, 5,25),c(30,25, 10,25))

for(n in 1:3){
      
  ################# Generating clusterised subsample ###################
  data1 <- runifpoint(numObs[[n]][1], win=W, giveup=1000, nsim=1, drop=TRUE, ex=NULL)
  
  data1.cords <- coords(data1)
  data1.cords.sp<-SpatialPoints(data1.cords)
  
  #setting the projection aligned with spatstat - plenar projection
  proj4string(data1.cords.sp)<-CRS("+proj=merc +datum=NAD83") 
  
  centrum1<- as.matrix(cbind(data1.cords.sp@coords[,1], data1.cords.sp@coords[,2]))
  
  # draw points around the centroid in a given radius 
  pointsMat1<-matrix(NA, nrow=0, ncol=2)
  for(i in 1:nrow(centrum1)){
    punkt<-centrum1[i,]
    ileClust<-numObs[[n]][2]
    
    rmaximum <- 10000
    angle <- runif(ileClust)*3.14*2
    radius <- runif(ileClust)*rmaximum
    
    xi <- punkt[1]+radius*cos(angle)
    yi <- punkt[2]+radius*sin(angle)
    
    pointsMat1<-rbind(pointsMat1, cbind(xi,yi))
    
  }
  
  pointsMat1.sp<-SpatialPoints(pointsMat1)
  #setting the projection aligned with spatstat - plenar projection
  proj4string(pointsMat1.sp)<-CRS("+proj=merc +datum=NAD83") 
  
  # check which points are within powiat Lublin 
  powLub <- inside.owin(pointsMat1.sp@coords[,1], pointsMat1.sp@coords[,2], W.lub)
  
  # change projection and create Spatial Points dataset
  pointsMat1.sp<-spTransform(pointsMat1.sp, CRS("+proj=longlat +datum=NAD83")) #sferyczne
  
  # creating a data frame with simulated companies characteristics
  data1<-cbind(pointsMat1.sp@coords, powLub)
  data1<-as.data.frame(data1)
  
  data1$roa_geo<-ifelse(data1$powLub==1, 1.5,0)
  data1$subsample<-4
  
  
  ################# Generating clusterised subsample v2 ###################
  data2 <- runifpoint(numObs[[n]][3], win=W, giveup=1000, nsim=1, drop=TRUE, ex=NULL)
  
  data2.cords <- coords(data2)
  data2.cords.sp<-SpatialPoints(data2.cords)
  
  #setting the projection aligned with spatstat - plenar projection
  proj4string(data2.cords.sp)<-CRS("+proj=merc +datum=NAD83") 
  
  centrum2<- as.matrix(cbind(data2.cords.sp@coords[,1], data2.cords.sp@coords[,2]))
  
  # draw points around the centroid in a given radius 
  pointsMat2<-matrix(NA, nrow=0, ncol=2)
  for(i in 1:nrow(centrum2)){
    punkt<-centrum2[i,]
    ileClust<-numObs[[n]][4]
    
    rmaximum <- 10000
    angle <- runif(ileClust)*3.14*2
    radius <- runif(ileClust)*rmaximum
    
    xi <- punkt[1]+radius*cos(angle)
    yi <- punkt[2]+radius*sin(angle)
    
    pointsMat2<-rbind(pointsMat2, cbind(xi,yi))
    
  }
  
  pointsMat2.sp<-SpatialPoints(pointsMat2)
  #setting the projection aligned with spatstat - plenar projection
  proj4string(pointsMat2.sp)<-CRS("+proj=merc +datum=NAD83") 
  
  # check which points are within powiat Lublin 
  powLub <- inside.owin(pointsMat2.sp@coords[,1], pointsMat2.sp@coords[,2], W.lub)
  
  # change projection and create Spatial Points dataset
  pointsMat2.sp<-spTransform(pointsMat2.sp, CRS("+proj=longlat +datum=NAD83")) #sferyczne
  
  # creating a data frame with simulated companies characteristics
  data2<-cbind(pointsMat2.sp@coords, powLub)
  data2<-as.data.frame(data2)
  
  data2$roa_geo<-ifelse(data2$powLub==1, 1.5,0)
  data2$subsample<-2
  
  #################################################################
  ############ Merging subsamples ############################
  
  fulldata<-rbind(data1, data2)
  fulldata<-as.data.frame(fulldata)
  fulldata$sector<-NA
  #fulldata[fulldata$subsample==1,]$sector<-'agri'
  fulldata[fulldata$subsample==2,]$sector<-'constr'
  fulldata[fulldata$subsample==4,]$sector<-'serv'
  
  fulldata$roa_sec<-NA
  #fulldata[fulldata$sector=='agri',]$roa_sec<-2
  #fulldata[fulldata$sector=='prod',]$roa_sec<-3.5
  fulldata[fulldata$sector=='constr',]$roa_sec<-5
  fulldata[fulldata$sector=='serv',]$roa_sec<-8
  
  
  library(spdep)
  # calculating euclinean distance between a point and Lublin's centroid
  coords<-as.matrix(data.frame(x=fulldata$x, y=fulldata$y)) 
  core<-c(22.5666700, 51.2500000) # Lublin's center coordinates
  fulldata$dist<-spDistsN1(coords, core, longlat=TRUE)
  
  # specifying roa_dist as a squared relation to the distance from the region's center
  fulldata$roa_dist<-0
  fulldata$roa_dist<- (1.5-(fulldata$dist/100))^2
  
  # adding all elements to create unique mean utilizable for generating the ROA value
  fulldata$roa_param<-fulldata$roa_sec+fulldata$roa_geo+fulldata$roa_dist 
  
  # generating ROA with individualized mean
  for(i in 1:nrow(fulldata)){
    fulldata$roa[i]<-rnorm(1, fulldata$roa_param[i], 0.045)}
  
  # binary values for sector recognition
  fulldata$agri<-ifelse(fulldata$sector=="agri",1,0)
  fulldata$prod<-ifelse(fulldata$sector=="prod",1,0)
  fulldata$constr<-ifelse(fulldata$sector=="constr",1,0)
  fulldata$serv<-ifelse(fulldata$sector=="serv",1,0)

  
  fulldata$sample <- numObs[[n]][1]*numObs[[n]][2]+numObs[[n]][3]*numObs[[n]][4]
  
  # calculating Clark-Evans statistics
  
  # creating a collection of coordinates
  fulldataBis<-SpatialPoints(cbind(fulldata$xi, fulldata$yi))
  proj4string(fulldataBis)<-CRS("+proj=longlat +datum=NAD83") 
  
  # changing projection to work with spatstat
  fulldataBis<-spTransform(fulldataBis, CRS("+proj=merc +datum=NAD83"))
  

  clarkx<-cbind(fulldataBis@coords)
  clarkx<-as.data.frame(clarkx)
  
  clarkx.ppp <- as.ppp(clarkx,W)
  
  clarkx.ppp <- as.ppp(clarkx, W)
  clarkResult <- clarkevans(clarkx.ppp)
  print(clarkResult)
  
  
  datax<-rbind(datax,fulldata)
  

}

datax2 <- datax[-1,]  

library(tidyverse)
library(sf)
library(osmdata)


datax.sf <- st_as_sf(datax2, coords = c("xi", "yi"), crs = "+proj=longlat +datum=NAD83")

woj.sf <- st_read("dane/wojewodztwa.shp")
woj.sf <- st_transform(woj.sf, crs = "+proj=longlat +datum=NAD83")

# limit to lubelskie
lub.woj.sf <- woj.sf %>% filter(jpt_nazwa_=='lubelskie') 

sampleLabs <- c("200 points", "500 points", "1000 points")
names(sampleLabs) <- c(200, 500, 1000)

ggplot() + 
  geom_sf(lub.woj.sf, mapping = aes(geometry = geometry))+
  geom_sf(datax.sf, mapping = aes(geometry = geometry, col = roa))+
  facet_grid(~sample, labeller = labeller(sample = sampleLabs))+
  #scale_colour_gradient(low = "#7b6afc", high ="#0b005e",  aesthetics = "colour")+
  scale_color_viridis_c()+
  labs(colour="ROA", title = 'Type F')+
  theme_minimal()+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))

# export for 950 width, 400 height