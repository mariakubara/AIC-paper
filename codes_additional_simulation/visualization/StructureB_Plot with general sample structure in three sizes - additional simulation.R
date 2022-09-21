library(tidyverse)

# Generating datasets

# setting the seed for reproducibility 
set.seed(601) 
#set.seed(12345678)

datax <- matrix(NA, ncol =15)
colnames(datax) <- c("x", "y", "powLub", "roa_geo", "sector", "roa_sec","dist", "roa_dist" ,
                     "roa_param", "roa", "agri", "prod", "constr", "serv", "sample" )

numObs<-c(200,500,1000)

for(n in 1:3){

  ################# Generating random points (centroids for clusters) ###################
  uniformSet <- runifpoint(numObs[n], win=W, giveup=1000, nsim=1, drop=TRUE, ex=NULL)
  
  data<-uniformSet
  data.cords <- coords(data)
  
  pointsMat.sp<-SpatialPoints(data.cords)
  proj4string(pointsMat.sp)<-CRS("+proj=merc +datum=NAD83") 
  plot(pointsMat.sp)
  
  ########## Check if points are withing the lubelskie window #####################
  
  # check which points are within powiat Lublin 
  # (geographical premium for being near the region's capital city)
  
  powLub <- inside.owin(pointsMat.sp@coords[,1], pointsMat.sp@coords[,2], W.lub)
  
  # change projection and create Spatial Points dataset
  pointsMat.sp<-spTransform(pointsMat.sp, CRS("+proj=longlat +datum=NAD83")) #sferyczne
  
  # creating a data frame with simulated companies characteristics
  data1<-cbind(pointsMat.sp@coords, powLub)
  data1<-as.data.frame(data1)
  
  # geographical premium
  data1$roa_geo<-ifelse(data1$powLub==1, 1.5,0)
  
  fulldata<-as.data.frame(data1)
  fulldata$sector<-'agri'
  
  fulldata$roa_sec<-NA
  fulldata[fulldata$sector=='agri',]$roa_sec<-2
  
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
  
  fulldata$sample <- numObs[n]
  
  # calculating Clark-Evans statistics
  
  # creating a collection of coordinates
  fulldataBis<-SpatialPoints(cbind(fulldata$x, fulldata$y))
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

datax.sf <- st_as_sf(datax2, coords = c("x", "y"), crs = "+proj=longlat +datum=NAD83")

sampleLabs <- c("200 points", "500 points", "1000 points")
names(sampleLabs) <- c(200, 500, 1000)

ggplot() + 
  geom_sf(lub.woj.sf, mapping = aes(geometry = geometry))+
  geom_sf(datax.sf, mapping = aes(geometry = geometry, col = roa))+
  facet_grid(~sample, labeller = labeller(sample = sampleLabs))+
  #scale_colour_gradient(low = "#7b6afc", high ="#0b005e",  aesthetics = "colour")+
  scale_color_viridis_c()+
  labs(colour="ROA", title = 'Type B')+
  theme_minimal()+
  theme(panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12))

# export for 950 width, 400 height