#loading our self-defined function code
source("codes_main_paper/semiVarianceKnn.R")

#setting the range of knn to be tested (2:maxKnn)
maxKnn <- 50

# preparing for the loop which will caluclate semivariance for the 
# initial simulation with clusters (with 20 observations each) in three sample sizes

numObs<-list(c(10,20),c(25,20),c(50,20))

resultsMat<-matrix(0, nrow=1, ncol=maxKnn) 

# vector with seeds that will be used to ensure replicability of a simulation
seed<-seq(1,1000,50)

for(n in 1:3){
  for(s in 1:20){ 
    # setting the seed for reproducibility 
    set.seed(seed[s])
    
    ### STEP 1 - generate the data structure (recognisable by our self-defined semi-variance function)
    
    ################# Generating random points (centroids for clusters) ###################
    uniformSet <- runifpoint(numObs[[n]][1], win=W, giveup=1000, nsim=1, drop=TRUE, ex=NULL)
    
    data<-uniformSet
    data.cords <- coords(data)
    data.cords.sp<-SpatialPoints(data.cords)
    
    #setting the projection aligned with spatstat - plenar projection
    proj4string(data.cords.sp)<-CRS("+proj=merc +datum=NAD83") 
    
    centroid<- as.matrix(cbind(data.cords.sp@coords[,1], data.cords.sp@coords[,2]))
    
    
    #2 draw points around the centroid in a given radius 
    pointsMat<-matrix(NA, nrow=0, ncol=2)
    
    for(i in 1:nrow(centroid)){
      point<-centroid[i,]
      howMany<-numObs[[n]][2]
      
      rmaximum <- 10000
      angle <- runif(howMany)*3.14*2
      radius <- runif(howMany)*rmaximum
      
      xi <- point[1]+radius*cos(angle)
      yi <- point[2]+radius*sin(angle)
      
      pointsMat<-rbind(pointsMat, cbind(xi,yi))
      
    }
    
    pointsMat.sp<-SpatialPoints(pointsMat)
    proj4string(pointsMat.sp)<-CRS("+proj=merc +datum=NAD83")
    
    ###############################    3    #########################################
    ########## Check if points are withing the lubelskie window #####################
    
    # check which points are within powiat Lublin 
    # (geographical premium for being near the region's capital city)
    
    powLub <- inside.owin(pointsMat.sp@coords[,1], pointsMat.sp@coords[,2], W.lub)
    
    # change projection and create Spatial Points dataset
    pointsMat.sp<-spTransform(pointsMat.sp, CRS("+proj=longlat +datum=NAD83")) 
    
    # creating a data frame with simulated companies characteristics
    data1<-cbind(pointsMat.sp@coords, powLub)
    data1<-as.data.frame(data1)
    
    # geographical premium
    data1$roa_geo<-ifelse(data1$powLub==1, 1.5,0)
    
    fulldata<-as.data.frame(data1)
    fulldata$sector<-'serv'
    
    fulldata$roa_sec<-NA
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
    
    datax<-fulldata
    
    crds<-as.matrix(cbind(datax$x, datax$y))
    
    ### STEP 2 - calculate semivariance for the data structure with specific 
    # coordinates, dependent variable and maximum Knn (range of knn 2:maxKnn)
    
    resul <- semiVarKnn(crds,datax$roa,maxKnn)
    
    # write results to the output file
    resultsMat<-rbind(resultsMat, resul)
    write.table(resultsMat, file="results_simulations/semiVarianceClusters2.csv", row.names = F, col.names = T)
    
  }
}